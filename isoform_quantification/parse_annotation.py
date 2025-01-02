#!/usr/bin/python
from operator import itemgetter, attrgetter
import re
import sys
import time
import numpy as np
import concurrent.futures
from util import sync_reference_name
import pickle
import copy
import config
from patch_mp import patch_mp_connection_bpo_17560
from parse_annotation_features import save_gene_structure_features
##########
def split_and_sort_exons(gene_exons_dict):
    new_gene_exons_dict = {}
    for chr_name in gene_exons_dict:
        new_gene_exons_dict[chr_name] = {}
        for gene_name in gene_exons_dict[chr_name]:
            new_exon_sorted = copy.deepcopy(gene_exons_dict[chr_name][gene_name])
            new_exon_sorted = sorted(new_exon_sorted,key=lambda x:(x[0],x[1]))
            i = 1
            while (i < len(new_exon_sorted)):
                if ((new_exon_sorted[i][0] == new_exon_sorted[i - 1][0]) and 
                    (new_exon_sorted[i][1] == new_exon_sorted[i - 1][1])):
                    del new_exon_sorted[i]  # Delete the repeated exon
                elif (new_exon_sorted[i][0] <= new_exon_sorted[i - 1][1]):
                    temp_exons = sorted([new_exon_sorted[i][0], new_exon_sorted[i - 1][0], new_exon_sorted[i][1], new_exon_sorted[i - 1][1]])
                    del new_exon_sorted[i - 1]  # Delete the two overlapping exons
                    del new_exon_sorted[i - 1]
                    if (temp_exons[0] == temp_exons[1]):  # Based on two exons overlap type, re-generate 2 or 3 new ones
                        new_exon_sorted.insert(i - 1, [temp_exons[1], temp_exons[2]])
                        new_exon_sorted.insert(i , [temp_exons[2] + 1, temp_exons[3]])
                    elif (temp_exons[2] == temp_exons[3]):
                        new_exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        new_exon_sorted.insert(i , [temp_exons[1], temp_exons[2]])
                    else:
                        new_exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        new_exon_sorted.insert(i, [temp_exons[1], temp_exons[2]])
                        new_exon_sorted.insert(i + 1, [temp_exons[2] + 1, temp_exons[3]])
                    new_exon_sorted = sorted(new_exon_sorted, key=itemgetter(0, 1))  # re-sort the exon positions
                else:
                    i += 1
            new_exon_sorted = sorted(new_exon_sorted, key=itemgetter(0, 1))
            # marginal_points = set()
            # for [start_pos,end_pos] in new_exon_sorted:
            #     if start_pos == end_pos:
            #         marginal_points.add(start_pos+1)
            # exon_sorted = []
            # for [start_pos,end_pos] in new_exon_sorted:
            #     if start_pos == end_pos:
            #         continue
            #     if start_pos in marginal_points:
            #         start_pos -= 1
            #     if end_pos in marginal_points:
            #         end_pos -= 1
            #     exon_sorted.append([start_pos,end_pos])
            # exon_sorted =  sorted(exon_sorted, key=itemgetter(0, 1))
            new_gene_exons_dict[chr_name][gene_name] = [[start_pos,end_pos,end_pos - start_pos + 1] for [start_pos,end_pos] in new_exon_sorted if start_pos < end_pos]
    return new_gene_exons_dict
def cal_region_length(region_name,point_coord_dict):
    assert '-' in region_name or ':' in region_name
    assert 'P' in list(point_coord_dict.keys())[0]
    if '-' not in region_name:
        exons = [region_name]
    else:
        exons = [region for region in region_name.split('-') if ':' in region]
    length = 0
    for exon in exons:
        points = exon.split(':')
        length += point_coord_dict[points[-1]] - point_coord_dict[points[0]] + 1
    return length
##########
def generate_exon_indicator_for_isoform_single_gene(args):
    single_gene_raw_isoform_exons_dict,exon_list,single_gene_gene_points_dict = args
    single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict = {},{}
    for isoform_name in single_gene_raw_isoform_exons_dict:
        isoform_exons = []
        #Question here?
#                 start_pos = [i-1 for i in single_gene_raw_isoform_exons_dict[isoform_name]['start_pos']]
        start_pos = single_gene_raw_isoform_exons_dict[isoform_name]['start_pos']
        end_pos = single_gene_raw_isoform_exons_dict[isoform_name]['end_pos']
        num_exons = len(single_gene_raw_isoform_exons_dict[isoform_name]['start_pos'])
        # Check the exon regions
        j = 0
        for i in range(len(exon_list)):
            flag = (j < num_exons)
            while (flag):
                p0 = exon_list[i][0]
                p1 = exon_list[i][1]
                l = exon_list[i][2]
                if ((p0 >= start_pos[j]) and
                    (p1 <= end_pos[j])):
                    region_name = 'P' + str(single_gene_gene_points_dict[p0]) + ':' + 'P' + str(single_gene_gene_points_dict[p1])
                    # if (l >= READ_LEN):  # A valid region for exon
                    temp_isoform_name = set()
                    temp_isoform_name.add(isoform_name)
                    if (region_name in single_gene_gene_regions_dict):
                        temp_isoform_name = temp_isoform_name.union(single_gene_gene_regions_dict[region_name])
                    single_gene_gene_regions_dict[region_name] = temp_isoform_name
                    isoform_exons.append(region_name)
                    flag = False
                elif (p0 <= start_pos[j]):
                    flag = False
                else:
                    j += 1
                    flag = (j < num_exons) 
        # Check the junction regions
        for i in range(len(exon_list) - 1):     # the last exon can not have any junction
            p0 = exon_list[i][0]
            p1 = exon_list[i][1]
            l  = exon_list[i][2]
            exon_region_name = 'P' + str(single_gene_gene_points_dict[p0]) + ':P' + str(single_gene_gene_points_dict[p1])
            if (exon_region_name not in isoform_exons):   # this exon is not in the isoform
                continue
            # if (l<READ_JUNC_MIN_MAP_LEN):
            #     continue
            region_name = 'P{}:P{}'.format(str(single_gene_gene_points_dict[p0]),str(single_gene_gene_points_dict[p1]))
            for j in range(i+1,len(exon_list)):
                p0_temp = exon_list[j][0]
                p1_temp = exon_list[j][1]
                l_temp  = exon_list[j][2]
                exon_region_name = 'P' + str(single_gene_gene_points_dict[p0_temp]) + ':P' + str(single_gene_gene_points_dict[p1_temp])
                if (exon_region_name not in isoform_exons):   # this exon is not in the isoform
                    continue
                # if (l_temp<READ_JUNC_MIN_MAP_LEN):
                #     continue
                if (region_name.split(':')[-1] == exon_region_name.split(':')[0]):
                    region_name += ':P{}'.format(str(single_gene_gene_points_dict[p1_temp]))
                else:
                    region_name += '-P{}:P{}'.format(str(single_gene_gene_points_dict[p0_temp]),str(single_gene_gene_points_dict[p1_temp]))
                temp_isoform_set = {isoform_name}
                if (region_name in single_gene_gene_regions_dict):
                    temp_isoform_set = temp_isoform_set.union(single_gene_gene_regions_dict[region_name])
                single_gene_gene_regions_dict[region_name] = temp_isoform_set
    single_gene_regions_len_dict = {}
    point_coord_dict = {}
    for p in single_gene_gene_points_dict:
        point_coord_dict['P{}'.format(single_gene_gene_points_dict[p])] = int(p)
    # filter out region with short ending exon
    # region_names = list(single_gene_gene_regions_dict.keys())
    # for region in region_names:
    #     if '-' in region:
    #         last_exon = region.split('-')[-1]
    #         last_exon_len = cal_region_length(last_exon,point_coord_dict)
    #         if last_exon_len < READ_JUNC_MIN_MAP_LEN:
    #             del single_gene_gene_regions_dict[region]
    for region in single_gene_gene_regions_dict:
        region_len = cal_region_length(region,point_coord_dict)
        single_gene_regions_len_dict[region] = region_len
        # if region is an exon
        if (region.count(':') == 1 and '-' not in region):
            for isoform_name in single_gene_gene_regions_dict[region]:
                if isoform_name in single_gene_isoforms_regions_len_dict:
                    single_gene_isoforms_regions_len_dict[isoform_name] += region_len
                else:
                    single_gene_isoforms_regions_len_dict[isoform_name] = region_len

    return single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict,single_gene_regions_len_dict
def generate_exon_indicator_for_isoform(gene_exons_dict,gene_points_dict,raw_isoform_exons_dict,threads):
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        isoforms_regions_len_dict[chr_name],gene_regions_dict[chr_name],genes_regions_len_dict[chr_name] = {},{},{}
        
    if threads == 1:
        for chr_name in raw_isoform_exons_dict:
            for gene_name in raw_isoform_exons_dict[chr_name]:
                isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = generate_exon_indicator_for_isoform_single_gene((raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name]))
    else:
        list_of_all_genes_chrs = [(gene_name,chr_name) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
        list_of_args = [(raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name]) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
        chunksize, extra = divmod(len(list_of_all_genes_chrs), threads)
        if extra:
            chunksize += 1
        print('[INFO] Using {} threads'.format(threads),flush=True)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, executor.map(generate_exon_indicator_for_isoform_single_gene,list_of_args,chunksize=chunksize)):
                isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = result
    # No parralize
    # for chr_name in raw_isoform_exons_dict:
    #     for gene_name in raw_isoform_exons_dict[chr_name]:
    #         isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = generate_exon_indicator_for_isoform_single_gene((raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name]))
    return isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict
#######################
def is_same_structure_isoform(raw_isoform_exons_dict,isoform_A,isoform_B):
    if len(raw_isoform_exons_dict[isoform_A]['start_pos']) == len(raw_isoform_exons_dict[isoform_B]['start_pos']):
        for i in range(len(raw_isoform_exons_dict[isoform_A]['start_pos'])):
            if raw_isoform_exons_dict[isoform_A]['start_pos'][i] != raw_isoform_exons_dict[isoform_B]['start_pos'][i] or raw_isoform_exons_dict[isoform_A]['end_pos'][i] != raw_isoform_exons_dict[isoform_B]['end_pos'][i]:
                return False
        return True
    else:
        return False
def parse_annotation(ref_annotation_path,threads):
    patch_mp_connection_bpo_17560()
    #gene_points_dict store the index of the point value in ascending order for each gene
    file_read = open(ref_annotation_path, 'r')
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = {},{},{},{},{}
    num_exons = 0
    for line in file_read:
        if line.lstrip()[0] == "#":
            continue
        fields = line.split('\t')
        if (fields[2] != 'exon'):
            continue
        chr_name = fields[0]
        converted_chr_name = sync_reference_name(fields[0])
        if (converted_chr_name.isnumeric()):
            chr_name = converted_chr_name
        # if chr_name == '':
        #     continue
        num_exons += 1
        gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
        isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
        # use 1 based index
        start_pos = int(fields[3])
        end_pos = int(fields[4])
        if start_pos > end_pos:
            continue
        # if chr_name in gene_exons_dict:
        #     if gene_name in gene_exons_dict[chr_name]:
        #         if (start_pos,end_pos) in gene_exons_dict[chr_name][gene_name]:
        #             continue
        #initialize dict
        if chr_name not in gene_exons_dict:
            gene_exons_dict[chr_name],gene_points_dict[chr_name],gene_isoforms_dict[chr_name],gene_isoforms_length_dict[chr_name],raw_isoform_exons_dict[chr_name] = {},{},{},{},{}
        if gene_name not in gene_exons_dict[chr_name]:
            gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],gene_isoforms_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],raw_isoform_exons_dict[chr_name][gene_name]= set(),{},[],{},{}
        gene_exons_dict[chr_name][gene_name].add((start_pos, end_pos))

        if isoform_name not in gene_isoforms_length_dict[chr_name][gene_name]:
            gene_isoforms_length_dict[chr_name][gene_name][isoform_name] = 0
            raw_isoform_exons_dict[chr_name][gene_name][isoform_name] = {'region_pos':[]}
            gene_isoforms_dict[chr_name][gene_name].append(isoform_name)
        
        # note here the base on both end included in our system
        gene_isoforms_length_dict[chr_name][gene_name][isoform_name] += end_pos - start_pos + 1
        raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos'].append([start_pos,end_pos])
    file_read.close()
    for chr_name in gene_exons_dict:
        for gene_name in gene_exons_dict[chr_name]:
            exon_list = []
            for (start_pos,end_pos) in gene_exons_dict[chr_name][gene_name]:
                exon_list.append([start_pos,end_pos])
            gene_exons_dict[chr_name][gene_name] = exon_list
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name].copy():
                region_pos = sorted(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos'],key=itemgetter(0, 1))
                raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'] = [start_pos for [start_pos,end_pos] in region_pos]
                raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'] = [end_pos for [start_pos,end_pos] in region_pos]
                del raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos']
    same_structure_isoform_dict = {}
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            isoform_names_to_join = set()
            for isoform_A in raw_isoform_exons_dict[chr_name][gene_name]:
                if isoform_A in isoform_names_to_join:
                    continue
                for isoform_B in raw_isoform_exons_dict[chr_name][gene_name]:
                    if (isoform_A != isoform_B) and (is_same_structure_isoform(raw_isoform_exons_dict[chr_name][gene_name],isoform_A,isoform_B)):
                        if config.same_struc_isoform_handling == 'merge':
                            isoform_names_to_join.add(isoform_A)
                            isoform_names_to_join.add(isoform_B)
                        # else:
                        #     if chr_name not in same_structure_isoform_dict:
                        #         same_structure_isoform_dict[chr_name] = {}
                        #     if gene_name not in same_structure_isoform_dict[chr_name]:
                        #         same_structure_isoform_dict[chr_name][gene_name] = {}
                        #     if isoform_A not in same_structure_isoform_dict[chr_name][gene_name]:
                        #         same_structure_isoform_dict[chr_name][gene_name][isoform_A] = set()

                        #     same_structure_isoform_dict[chr_name][gene_name][isoform_A].add(isoform_B)
                        #     isoform_names_to_join.add(isoform_B)
            if config.same_struc_isoform_handling == 'merge':
                if len(isoform_names_to_join) > 0:
                    joined_isoform_name = '-'.join(list(isoform_names_to_join))
                    original_isoform_name = list(isoform_names_to_join)[0]
                    raw_isoform_exons_dict[chr_name][gene_name][joined_isoform_name] = copy.deepcopy(raw_isoform_exons_dict[chr_name][gene_name][original_isoform_name])
                    gene_isoforms_length_dict[chr_name][gene_name][joined_isoform_name]= copy.deepcopy(gene_isoforms_length_dict[chr_name][gene_name][original_isoform_name])
                    gene_isoforms_dict[chr_name][gene_name].append(joined_isoform_name)
                for isoform_name in isoform_names_to_join:
                    del raw_isoform_exons_dict[chr_name][gene_name][isoform_name]
                    del gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                    gene_isoforms_dict[chr_name][gene_name].remove(isoform_name)
    raw_gene_exons_dict = gene_exons_dict.copy()
    for chr_name in gene_isoforms_dict:
        for gene_name in list(gene_isoforms_dict[chr_name].keys()):
            if len(gene_isoforms_dict[chr_name][gene_name]) == 0 :
                del gene_exons_dict[chr_name][gene_name]
                del raw_isoform_exons_dict[chr_name][gene_name]
                del gene_isoforms_length_dict[chr_name][gene_name]
                del gene_isoforms_dict[chr_name][gene_name]
                del gene_points_dict[chr_name][gene_name]
    raw_gene_exons_dict = copy.deepcopy(gene_exons_dict)
    gene_exons_dict = split_and_sort_exons(gene_exons_dict)
    # index the point position
    for chr_name in gene_exons_dict:
        for gene_name in gene_exons_dict[chr_name]:
            point_index = 0
            for [start_pos,end_pos,_] in gene_exons_dict[chr_name][gene_name]:
                for pos in [start_pos,end_pos]:
                    if pos not in gene_points_dict[chr_name][gene_name]:
                        gene_points_dict[chr_name][gene_name][pos] = point_index
                        point_index += 1
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = generate_exon_indicator_for_isoform(gene_exons_dict, gene_points_dict, raw_isoform_exons_dict,threads)
    save_gene_structure_features(gene_exons_dict,raw_gene_exons_dict,gene_isoforms_length_dict,genes_regions_len_dict,config.output_path)
    return [gene_exons_dict,gene_points_dict, gene_isoforms_dict,genes_regions_len_dict,
            isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict,same_structure_isoform_dict]
