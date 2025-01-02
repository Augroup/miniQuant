from parse_alignment import map_read,parse_read_line
from patch_mp import patch_mp_connection_bpo_17560
from parse_annotation_main import check_valid_region
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import numpy as np
import time
import random
import multiprocessing as mp
import os
from util import check_region_type
import config
import pickle
from pathlib import Path
import pandas as pd
import shutil
import glob

# from memory_profiler import profile
# def parse_alignment_iteration(alignment_file_path,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
#                     start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#                     READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,map_f,line_nums):
def debuginfoStr(info):
    print(info,flush=True)
    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]
    mem = int(memusage.strip())/1024
    print('Mem consumption: '+str(mem),flush=True)
def get_reads_isoform_info(output_path,gene_regions_read_mapping,LR_gene_regions_dict,LR_feature_dict):
    with open(f'{output_path}/temp/LR_alignments_dict/isoform_dict','rb') as f:
        [isoform_len_dict,isoform_exon_dict,strand_dict] = pickle.load(f)
    reads_isoform_info = {}
    expression_dict = {}
    unique_mapping_expression_dict = {}
    for rname in gene_regions_read_mapping:
        for gname in gene_regions_read_mapping[rname]:
            gene_unique_map_read_length = []
            gene_multi_map_read_length = []
            gene_3_end_truncation = []
            gene_5_end_truncation = []
            for region in gene_regions_read_mapping[rname][gname]:
                for read_mapping in gene_regions_read_mapping[rname][gname][region]:
                    read_start,read_end = read_mapping['read_pos']
                    read_name = read_mapping['read_name']
                    if read_name not in reads_isoform_info:
                        reads_isoform_info[read_name] = {}
                    for isoform in LR_gene_regions_dict[rname][gname][region]: 
                        # assign non-zero expression isoform
                        if isoform not in expression_dict:
                            expression_dict[isoform] = 0
                            unique_mapping_expression_dict[isoform] = 0
                        isoform_len = isoform_len_dict[isoform]
                        start_offset = 0
                        for [exon_start,exon_end] in isoform_exon_dict[isoform]:
                            if exon_end > read_start:
                                start_offset += read_start - exon_start + 1
                                break
                            else:
                                start_offset += exon_end - exon_start + 1
                        end_offset = isoform_len_dict[isoform] - start_offset - read_mapping['read_length']
                        if start_offset < 0:
                            start_offset = 0
                        if end_offset < 0:
                            end_offset = 0
                        if strand_dict[gname]  == '+':
                            reads_isoform_info[read_name][isoform] = [isoform_len,start_offset,end_offset]
                        else:
                            reads_isoform_info[read_name][isoform] = [isoform_len,end_offset,start_offset]
                    if len(LR_gene_regions_dict[rname][gname][region]) == 1:
                         # use unique mapping reads as initial expression
                        isoform = next(iter(LR_gene_regions_dict[rname][gname][region]))
                        expression_dict[isoform] += 1
                        unique_mapping_expression_dict[isoform] += 1
                        gene_unique_map_read_length.append(read_mapping['read_length'])
                        if strand_dict[gname]  == '+':
                            gene_5_end_truncation.append(start_offset/isoform_len)
                            gene_3_end_truncation.append(end_offset/isoform_len)
                        elif strand_dict[gname]  == '-':
                            gene_3_end_truncation.append(start_offset/isoform_len)
                            gene_5_end_truncation.append(end_offset/isoform_len)
                        else:
                            raise Exception('Stand wrong!')
                    else:
                        gene_multi_map_read_length.append(read_mapping['read_length'])
                        num_isoforms = len(LR_gene_regions_dict[rname][gname][region])
                        if config.inital_theta not in ['SR_unique','LR_unique','hybrid_unique']:
                            for isoform in LR_gene_regions_dict[rname][gname][region]:
                                expression_dict[isoform] += 1/num_isoforms
            if gname not in LR_feature_dict:
                LR_feature_dict[gname] = {'unique_mapping':gene_unique_map_read_length,'multi_mapping':gene_multi_map_read_length,'3_end_truncation':gene_3_end_truncation,'5_end_truncation':gene_5_end_truncation}
            else:
                LR_feature_dict[gname]['unique_mapping'] += gene_unique_map_read_length
                LR_feature_dict[gname]['multi_mapping'] += gene_multi_map_read_length
                LR_feature_dict[gname]['3_end_truncation'] += gene_3_end_truncation
                LR_feature_dict[gname]['5_end_truncation'] += gene_5_end_truncation
    # print(len(reads_isoform_info))
    return reads_isoform_info,expression_dict,unique_mapping_expression_dict,LR_feature_dict
def get_read_len_dist(reads_isoform_info):
    read_len_dict = {}
    read_len_dist_dict = {}
    for read in reads_isoform_info:
        max_read_len  = 0
        for isoform,[isoform_len,f_end,t_end] in reads_isoform_info[read].items():
            read_len = isoform_len -  f_end - t_end
            if read_len > max_read_len:
                max_read_len = read_len
        read_len_dict[read] = max_read_len
        if max_read_len not in read_len_dist_dict:
            read_len_dist_dict[max_read_len] = 0
        read_len_dist_dict[max_read_len] += 1
    # read_len_dist = pd.Series(read_len_dist_dict).to_frame()
    # read_len_dist.columns = ['PDF']
    # read_len_dist = read_len_dist.sort_index(ascending=True)
    return read_len_dict,read_len_dist_dict
def parse_alignment_iteration(alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,CHR_LIST,output_path,aln_line_marker):
    os.nice(10)
    Path(f'{output_path}/temp/LR_alignments/').mkdir(exist_ok=True,parents=True)
    start_file_pos,end_file_pos,worker_id = aln_line_marker
    batch_id = 0
    points_dict,interval_tree_dict, gene_regions_dict,\
        start_pos_list, start_gname_list, end_pos_list, end_gname_list = {},{},{},{},{},{},{}
    LR_feature_dict = {}
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_pos = {}
        local_gene_regions_read_count = {}
        local_gene_regions_read_length = {}
        aln_file.seek(start_file_pos)
        max_buffer_size = 1e5
        buffer_size = 0
        while True:
            line = aln_file.readline()
            # if aln_file.tell() > end_file_pos or not line:
            #     break
            if aln_file.tell() > end_file_pos or not line:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'):
                    continue
                aln_line = parse_read_line(line)
                [_, _, rname, _] = aln_line
                if rname not in points_dict:
                    points_dict,interval_tree_dict,start_pos_list, start_gname_list, end_pos_list, end_gname_list = {},{},{},{},{},{}
                    with open(f'{output_path}/temp/LR_alignments_dict/{rname}','rb') as f:
                        [points_dict[rname],interval_tree_dict[rname],gene_regions_dict[rname],start_pos_list[rname],end_pos_list[rname],start_gname_list[rname],end_gname_list[rname]] = pickle.load(f)
                mapping = map_f(points_dict,interval_tree_dict, gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
                if (mapping['read_mapped']):
                    #random.seed(mapping['read_name'])
                    for mapping_area in mapping['mapping_area']:
                    # for mapping_are in [random.choice(mapping['mapping_area'])]:
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']

                        if rname not in local_gene_regions_read_count:
                            local_gene_regions_read_pos[rname] = {}
                            local_gene_regions_read_count[rname],local_gene_regions_read_length[rname] = {},{}
                        if gname not in local_gene_regions_read_count[rname]:
                            local_gene_regions_read_pos[rname][gname] = {}
                            local_gene_regions_read_count[rname][gname],local_gene_regions_read_length[rname][gname] = {},{}
                        if region_name not in local_gene_regions_read_count[rname][gname]:
                            local_gene_regions_read_pos[rname][gname][region_name] = []
                            local_gene_regions_read_count[rname][gname][region_name],local_gene_regions_read_length[rname][gname][region_name] = 0,[]
                        local_gene_regions_read_pos[rname][gname][region_name].append(mapping)
                        local_gene_regions_read_count[rname][gname][region_name] += 1 
                        local_gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
                        buffer_size += 1
            except Exception as e:
                # tb = traceback.format_exc()
                # print(Exception('Failed to on ' + line, tb))
                continue
            if buffer_size > max_buffer_size:
                reads_isoform_info,expression_dict,unique_mapping_expression_dict,LR_feature_dict = get_reads_isoform_info(output_path,local_gene_regions_read_pos,gene_regions_dict,LR_feature_dict)
                read_len_dict,read_len_dist_dict = get_read_len_dist(reads_isoform_info)
                with open(f'{output_path}/temp/LR_alignments/reads_{worker_id}_{batch_id}','wb') as f:
                    pickle.dump([reads_isoform_info,read_len_dict],f)
                with open(f'{output_path}/temp/LR_alignments/dist_{worker_id}_{batch_id}','wb') as f:
                    pickle.dump([read_len_dist_dict,expression_dict],f)
                with open(f'{output_path}/temp/LR_alignments/read_count_length_{worker_id}_{batch_id}','wb') as f:
                    pickle.dump([local_gene_regions_read_count,local_gene_regions_read_length],f)
                with open(f'{output_path}/temp/machine_learning/LR_feature_dict_{worker_id}_{batch_id}','wb') as f:
                    pickle.dump(LR_feature_dict,f)
                batch_id += 1
                del reads_isoform_info
                del expression_dict
                del unique_mapping_expression_dict
                local_gene_regions_read_pos = {}
                local_gene_regions_read_count = {}
                local_gene_regions_read_length = {}
                LR_feature_dict = {}
                buffer_size = 0
        if buffer_size > 0:
            reads_isoform_info,expression_dict,unique_mapping_expression_dict,LR_feature_dict =get_reads_isoform_info(output_path,local_gene_regions_read_pos,gene_regions_dict,LR_feature_dict)
            read_len_dict,read_len_dist_dict = get_read_len_dist(reads_isoform_info)
            with open(f'{output_path}/temp/LR_alignments/reads_{worker_id}_{batch_id}','wb') as f:
                pickle.dump([reads_isoform_info,read_len_dict],f)
            with open(f'{output_path}/temp/LR_alignments/dist_{worker_id}_{batch_id}','wb') as f:
                pickle.dump([read_len_dist_dict,expression_dict],f)
            with open(f'{output_path}/temp/LR_alignments/read_count_length_{worker_id}_{batch_id}','wb') as f:
                pickle.dump([local_gene_regions_read_count,local_gene_regions_read_length],f)
            with open(f'{output_path}/temp/machine_learning/LR_feature_dict_{worker_id}_{batch_id}','wb') as f:
                pickle.dump(LR_feature_dict,f)
            local_gene_regions_read_count = {}
            local_gene_regions_read_length = {}
            LR_feature_dict = {}
    return 
# @profile
# def get_aln_line_marker(alignment_file_path,threads):
#     with open(alignment_file_path, 'r') as aln_file:
#         line_offset = []
#         offset = 0
#         for line in aln_file:
#             if line[0] != '@':
#                 line_offset.append(offset)
#             offset += len(line)
#     num_aln_lines = len(line_offset)
#     chunksize, extra = divmod(num_aln_lines, threads)
#     if extra:
#         chunksize += 1
#     aln_line_marker = []
#     for i in range(threads):
#         aln_line_marker.append((line_offset[i*chunksize],chunksize,i))
#     return aln_line_marker
def get_aln_line_marker(alignment_file_path,threads):
    '''
    Split the sam file into THREADS chunks
    !Split by read
    '''
    with open(alignment_file_path,'r') as f:
        while True:
            line = f.readline()
            if line[0] != '@':
                header_size = f.tell()
                break
    file_stats = os.stat(alignment_file_path)
    total_bytes = file_stats.st_size-header_size
    chunksize, extra = divmod(total_bytes, threads)
    if extra:
        chunksize += 1
    byte_marker = []
    with open(alignment_file_path,'r') as f:
        for i in range(threads):
            if i == 0:
                start_offset = 0
                for line in f:
                    if line[0] != '@':
                        break
                    start_offset += len(line)
            else:
                f.seek(header_size+i*chunksize)
                f.readline()
                start_offset = f.tell()
            byte_marker.append(start_offset)
    byte_marker.append(total_bytes+header_size)
    aln_line_marker = []
    for i in range(threads):
         aln_line_marker.append((byte_marker[i],byte_marker[i+1],i))
    return aln_line_marker
def parse_alignment_EM(alignment_file_path,READ_JUNC_MIN_MAP_LEN,output_path,threads,CHR_LIST):
    # print(threads)
    patch_mp_connection_bpo_17560()
    # Create sorted end and start positions
    map_f = map_read
    # debuginfoStr('Before MP mem usage')
    pool = mp.Pool(threads)
    # partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path)
    # temp_queue = manager.Queue()    
    # watcher = pool.apply_async(mapping_listener, args=(temp_queue,gene_regions_read_count,gene_regions_read_length,gene_regions_read_pos))
    partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,CHR_LIST,output_path)
    futures = []
    aln_line_marker = get_aln_line_marker(alignment_file_path,threads)
    for marker in aln_line_marker:
        futures.append(pool.apply_async(partial_read_alignment,(marker,)))
    # with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #     futures = [executor.submit(partial_read_alignment,marker) for marker in aln_line_marker]
    #     concurrent.futures.wait(futures)
    for future in futures:
        future.get()
    # temp_queue.put('kill')
    # gene_regions_read_count,gene_regions_read_length,num_mapped_lines,gene_regions_read_pos = watcher.get()
    pool.close()
    pool.join()
def get_region_read_count_length(gene_regions_dict,output_path):
    gene_regions_read_count = {}
    gene_regions_read_length = {}
    for fpath in glob.glob(f'{output_path}/temp/LR_alignments/read_count_length_*_*'):
        with open(fpath,'rb') as f:
            [local_gene_regions_read_count,local_gene_regions_read_length] = pickle.load(f)
            for rname in local_gene_regions_read_count:
                if rname not in gene_regions_read_count:
                    gene_regions_read_count[rname] = {}
                    gene_regions_read_length[rname] = {}
                for gname in local_gene_regions_read_count[rname]:
                    if gname not in gene_regions_read_count[rname]:
                        gene_regions_read_count[rname][gname] = {}
                        gene_regions_read_length[rname][gname] = {}
                    for region in local_gene_regions_read_count[rname][gname]:
                        if region not in gene_regions_read_count[rname][gname]:
                            gene_regions_read_count[rname][gname][region] = 0
                            gene_regions_read_length[rname][gname][region] = []
                        gene_regions_read_count[rname][gname][region] += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_length[rname][gname][region] += local_gene_regions_read_length[rname][gname][region]
                
    read_lens = []
    for rname in gene_regions_read_length:
        for gname in gene_regions_read_length[rname]:
            for region in gene_regions_read_length[rname][gname]:
                read_lens += gene_regions_read_length[rname][gname][region]
        
    num_long_reads = 0
    gene_full_length_region_dict = defaultdict(lambda:{})
    isoform_max_num_exons_dict = {}
    for rname in gene_regions_dict:
        for gname in gene_regions_dict[rname]:
            regions_set = set()
            isoform_region_dict = defaultdict(lambda:set())
            for region in gene_regions_dict[rname][gname]:
                for isoform in gene_regions_dict[rname][gname][region]:
                    isoform_region_dict[isoform].add(region)
            for isoform in isoform_region_dict:
                max_region_exon_num = 0
                longest_region = ''
                for region in isoform_region_dict[isoform]:
                    region_exon_num = region.count(':')
                    if max_region_exon_num < region_exon_num:
                        max_region_exon_num = region_exon_num
                        longest_region = region
                isoform_max_num_exons_dict[isoform] = max_region_exon_num
                for region in isoform_region_dict[isoform]:
                    region_exon_num = region.count(':')
                    if region_exon_num == max_region_exon_num:
                        regions_set.add(region)
            gene_full_length_region_dict[rname][gname] = regions_set
    filtered_gene_regions_read_length = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
    for rname in gene_regions_read_length.copy():
        for gname in gene_regions_read_length[rname].copy():    
            # region_lens = []
            for region in gene_regions_read_length[rname][gname].copy():
                region_exon_num = region.count(':')
                read_length_list = []
                for read_length in gene_regions_read_length[rname][gname][region]:
                    is_valid_read = True
                    if (is_valid_read):
                        read_length_list.append(read_length)
                gene_regions_read_length[rname][gname][region] = read_length_list
                gene_regions_read_count[rname][gname][region] = len(read_length_list)
                num_long_reads += len(read_length_list)

                if gene_regions_read_count[rname][gname][region] == 0:
                    if region not in gene_full_length_region_dict[rname][gname]:
                        del gene_regions_read_length[rname][gname][region]
                        del gene_regions_read_count[rname][gname][region]
    return gene_regions_read_count,gene_regions_read_length,sum(read_lens),num_long_reads
   