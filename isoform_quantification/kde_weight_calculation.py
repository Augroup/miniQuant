
import sys
import numpy as np
import copy
from pathlib import Path
import pickle
import concurrent.futures
import os, shutil
import time
# import joblib
# import config
# model_path = '/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/human_simulation/data/model/human_NA12878_dRNA_Bham1_guppy/training_aligned_region_2d.pkl'
# kde_aligned = joblib.load(model_path)
def check_num_exons(region):
    return region.count(':')
def calculate_min_region_length(region_len_dict,READ_JUNC_MIN_MAP_LEN):
    region_min_length_dict = {}
    for region_name in region_len_dict:
        if '-' not in region_name:
            region_len = cal_exon_region_min_len(region_name,region_len_dict,READ_JUNC_MIN_MAP_LEN)
        else:
            exons = region_name.split('-')
            first_exon = exons[0]
            end_exon = exons[-1]
            if len(exons) > 2:
                inner_exons = exons[1:-1]
                inner_region_len = region_len_dict['-'.join(inner_exons)]
            else:
                inner_region_len = 0
            region_len = (region_len_dict[first_exon],region_len_dict[end_exon],inner_region_len)
#             first_exon_min_map_len = cal_exon_region_min_len(first_exon,region_len_dict,READ_JUNC_MIN_MAP_LEN,'end')
#             if first_exon_min_map_len > region_len_dict[first_exon]:
#                 first_exon_min_map_len = region_len_dict[first_exon]
#             end_exon_min_map_len = cal_exon_region_min_len(first_exon,region_len_dict,READ_JUNC_MIN_MAP_LEN,'start')
#             if end_exon_min_map_len > region_len_dict[end_exon]:
#                 end_exon_min_map_len = region_len_dict[end_exon]
#             region_len = first_exon_min_map_len + inner_region_len + end_exon_min_map_len,
        region_min_length_dict[region_name] = region_len
    return region_min_length_dict
def cal_exon_region_min_len(region_name,region_len_dict,READ_JUNC_MIN_MAP_LEN,junc_end=None):
    assert '-' not in region_name
    if check_num_exons(region_name) == 1:
        return READ_JUNC_MIN_MAP_LEN
    if check_num_exons(region_name) == 2:
        return READ_JUNC_MIN_MAP_LEN
    if check_num_exons(region_name) > 2:
        inner_region = ':'.join(region_name.split(':')[1:-1])
        inner_region_len = region_len_dict[inner_region]
        if junc_end is None:
            first_exon_len = 1
            end_exon_len = 1
            return inner_region_len +first_exon_len +  end_exon_len
        elif junc_end == 'end':
            first_exon_len = 1
            end_exon = ':'.join(region_name.split(':')[-2:])
            end_exon_len = region_len_dict[end_exon]
            return inner_region_len +first_exon_len +  end_exon_len
        elif junc_end == 'start':
            end_exon_len = 1
            first_exon = ':'.join(region_name.split(':')[:2])
            first_exon_len = region_len_dict[first_exon]
            return inner_region_len +first_exon_len +  end_exon_len
def cal_isoform_region_weight_single_gene(args):
    # isoform_region_weight_dict = {}
    # for rname in LR_gene_regions_dict:
    #     isoform_region_weight_dict[rname] = {}
    #     for gname in LR_gene_regions_dict[rname]:
    #         isoform_region_weight_dict[rname][gname] = {}
    [output_dir,rname,gname,READ_JUNC_MIN_MAP_LEN] = args
    READ_JUNC_MIN_MAP_LEN = 1
    start_time = time.time()
    with open(f'{output_dir}/temp/{rname}_{gname}.pkl','rb') as f:
        [LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict] = pickle.load(f)
    isoform_region_weight_dict = {}

    min_region_len_dict = calculate_min_region_length(LR_genes_regions_len_dict,READ_JUNC_MIN_MAP_LEN)
    num_regions = len(LR_gene_regions_dict)
    for region in LR_gene_regions_dict:
        max_region_len = LR_genes_regions_len_dict[region]
        min_region_len = min_region_len_dict[region]
        if '-' in region:
            region_type = 'junction'
        else:
            region_type = 'exon'
            if max_region_len < READ_JUNC_MIN_MAP_LEN:
                min_region_len = max_region_len
#             assert min_region_len <= max_region_len
#         xrange = range(min_region_len,max_region_len+1)
        for isoform in LR_gene_regions_dict[region]:
            if isoform not in isoform_region_weight_dict:
                isoform_region_weight_dict[isoform] = {}
            isoform_length = gene_isoforms_length_dict[isoform]
            if region_type == 'exon':
                prob_sum = 0
                xrange = range(min_region_len,max_region_len+1)
                X = np.array([[isoform_length,i] for i in xrange])
                X_prob = np.exp(kde_aligned.score_samples(X))
                for prob,t in zip(X_prob,range(min_region_len,max_region_len+1)):
                    prob_sum += prob * (max_region_len - t + 1)
            else:
                prob_sum = 0
                first_exon_len,end_exon_len,inner_region_len = min_region_len
                lb = inner_region_len+max(first_exon_len,end_exon_len)+READ_JUNC_MIN_MAP_LEN + 1
                ub = first_exon_len+end_exon_len+inner_region_len
                xrange = range(lb,ub+1)
                X = np.array([[isoform_length,i] for i in xrange])
                X_prob = np.exp(kde_aligned.score_samples(X))
                for prob,t in zip(X_prob,range(lb,ub+1)):
                    prob_sum += prob * (first_exon_len+end_exon_len+inner_region_len - t + 1)
                lb = inner_region_len+min(first_exon_len,end_exon_len)+READ_JUNC_MIN_MAP_LEN
                ub = inner_region_len+max(first_exon_len,end_exon_len)+READ_JUNC_MIN_MAP_LEN
                xrange = range(lb,ub+1)
                X = np.array([[isoform_length,i] for i in xrange])
                X_prob = np.exp(kde_aligned.score_samples(X))
                for prob,t in zip(X_prob,range(lb,ub+1)):
                    prob_sum += prob * min(first_exon_len,end_exon_len)
                lb = inner_region_len + 2 * READ_JUNC_MIN_MAP_LEN
                ub = inner_region_len+min(first_exon_len,end_exon_len)+READ_JUNC_MIN_MAP_LEN - 1
                xrange = range(lb,ub+1)
                X = np.array([[isoform_length,i] for i in xrange])
                X_prob = np.exp(kde_aligned.score_samples(X))
                for prob,t in zip(X_prob,range(lb,ub+1)):
                    prob_sum += prob * (t - 1 - inner_region_len)
            isoform_region_weight_dict[isoform][region] = prob_sum
    return isoform_region_weight_dict
def dump_to_disk(output_dir,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict):
    for rname in LR_gene_regions_dict:
        for gname in LR_gene_regions_dict[rname]:
            all_dict = [LR_gene_regions_dict[rname][gname],LR_genes_regions_len_dict[rname][gname],gene_isoforms_length_dict[rname][gname]]
            Path(f'{output_dir}/temp/').mkdir(exist_ok=True,parents=True)
            with open(f'{output_dir}/temp/{rname}_{gname}.pkl','wb') as f:
                pickle.dump(all_dict,f)
def multi_thread_management(list_of_all_genes_chrs,READ_JUNC_MIN_MAP_LEN,output_dir,threads):
    isoform_region_weight_dict = {}
    list_of_args = [[output_dir,chr_name,gene_name,READ_JUNC_MIN_MAP_LEN] for gene_name,chr_name in list_of_all_genes_chrs]
    chunksize, extra = divmod(len(list_of_args), threads)
    if extra:
        chunksize += 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, executor.map(cal_isoform_region_weight_single_gene,list_of_args,chunksize=chunksize)):
            if chr_name not in isoform_region_weight_dict:
                isoform_region_weight_dict[chr_name] = {}
            isoform_region_weight_dict[chr_name][gene_name] = result
    return isoform_region_weight_dict
def cal_isoform_region_weight(gene_regions_dict,gene_region_len_dict,gene_isoforms_length_dict,READ_JUNC_MIN_MAP_LEN,output_dir,threads):
    dump_to_disk(output_dir,gene_regions_dict,gene_region_len_dict,gene_isoforms_length_dict)
    list_of_all_genes_chrs = []
    for chr_name in gene_regions_dict:
        for gene_name in gene_regions_dict[chr_name]:
            # if gene_name in config.gene_list:
            list_of_all_genes_chrs.append((gene_name,chr_name))
    isoform_region_weight_dict =  multi_thread_management(list_of_all_genes_chrs,READ_JUNC_MIN_MAP_LEN,output_dir,threads)
    temp_folder = f'{output_dir}/temp/'
    for filename in os.listdir(temp_folder):
        file_path = os.path.join(temp_folder, filename)
        if '.pkl' not in file_path:
            continue
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('[INFO] Failed to delete %s. Reason: %s' % (file_path, e))
    return isoform_region_weight_dict

def create_new_matrix_dict(long_read_gene_matrix_dict,long_reads_isoform_region_weight_matrix_dict):
    new_long_read_gene_matrix_dict = long_read_gene_matrix_dict.copy()
    isoform_region_matrix = new_long_read_gene_matrix_dict['isoform_region_matrix']
    isoform_region_matrix[isoform_region_matrix != 0] = -1
    for isoform_name in long_read_gene_matrix_dict['isoform_names_indics']:
        if isoform_name not in long_reads_isoform_region_weight_matrix_dict:
            for split_isoform in isoform_name.split('-'):
                if split_isoform in long_reads_isoform_region_weight_matrix_dict:
                    isoform_index = long_read_gene_matrix_dict['isoform_names_indics'][isoform_name]
                    isoform_name = split_isoform
                    break
        else:
            isoform_index = long_read_gene_matrix_dict['isoform_names_indics'][isoform_name]
        if isoform_name in long_reads_isoform_region_weight_matrix_dict:
            for region in long_reads_isoform_region_weight_matrix_dict[isoform_name]:
                if region not in long_read_gene_matrix_dict['region_names_indics']:
                    continue
                region_index = long_read_gene_matrix_dict['region_names_indics'][region]
                if isoform_region_matrix[region_index][isoform_index] == -1:
                    isoform_region_matrix[region_index][isoform_index] = long_reads_isoform_region_weight_matrix_dict[isoform_name][region]
    isoform_region_matrix[isoform_region_matrix == -1] = 0
    sum_A = isoform_region_matrix.sum(axis=0)
    sum_A[sum_A==0] = 1
    isoform_region_matrix = isoform_region_matrix/sum_A
    new_long_read_gene_matrix_dict['isoform_region_matrix'] = isoform_region_matrix
    return new_long_read_gene_matrix_dict     