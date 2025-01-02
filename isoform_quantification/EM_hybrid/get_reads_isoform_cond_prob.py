import pandas as pd
import pickle
from pathlib import Path
import multiprocessing as mp
import numpy as np
import glob
import scipy
import scipy.sparse
import config
from EM_hybrid.util import convert_dict_to_sparse_matrix
def get_cdf(read_len_dist,l):
    return read_len_dist.loc[read_len_dist.index<l,'PDF'].sum()/read_len_dist['PDF'].sum()
def get_Sm_dict(read_len_dist,isoform_df):
    descending_read_len_dist = read_len_dist.sort_index(ascending=False)
    Sm_dict = {}
    previous_Sm= 0
    previous_read_len = 0
    first_read_len = True
    for read_len in descending_read_len_dist.index:
        if first_read_len:
            sub_isoform_df = isoform_df[(isoform_df['isoform_len']>=read_len)].copy()
        else:
            sub_isoform_df = isoform_df[(isoform_df['isoform_len']>=read_len) & (isoform_df['isoform_len']<previous_read_len)].copy()
        sub_isoform_len_set = set(sub_isoform_df['isoform_len'].values)
        sub_isoform_len_cdf_dict = {}
        for isoform_len in sub_isoform_len_set:
            sub_isoform_len_cdf_dict[isoform_len] = get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1) 
        if sub_isoform_df.shape[0] != 0:
            sub_isoform_df['sum'] = sub_isoform_df.apply(lambda rows:rows['theta']/sub_isoform_len_cdf_dict[rows['isoform_len']],axis=1)
            Lm = sub_isoform_df['sum'].sum()
        else:
            Lm = 0
        if first_read_len:
            Sm_dict[read_len] = Lm
            first_read_len = False
        else:
            Sm_dict[read_len] = previous_Sm + Lm
        previous_Sm = Sm_dict[read_len]
        previous_read_len = read_len
    return Sm_dict
def get_all_reads_isoform_cond_prob_LIQA_modified(args):
    worker_id,output_path,isoform_df,isoform_index_dict,read_len_dist,Sm_dict = args
    MIN_PROB = 1e-100
    num_batches = 0
    all_cond_prob_matrix = []
    for fpath in sorted(glob.glob(f'{output_path}/temp/LR_alignments/reads_{worker_id}_*')):
        batch_id = fpath.split('/')[-1].split('_')[-1]
        with open(fpath,'rb') as f:
            [reads_isoform_info,read_len_dict] = pickle.load(f)
        all_reads_isoform_cond_prob = {}
        read_index = -1
        for read in reads_isoform_info:
            if len(reads_isoform_info[read]) != 0:
                read_index += 1
            read_len = read_len_dict[read]
            for isoform in reads_isoform_info[read]:
                isoform_len = isoform_df.loc[isoform,'isoform_len']
                pdf = read_len_dist.loc[read_len,'PDF']/read_len_dist['PDF'].sum()
                cdf = get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1)
                Sm = Sm_dict[read_len]
                if cdf == 0 or Sm == 0:
                    cond_prob = pdf
                else:
                    if config.LR_cond_prob_calc == 'form_1':
                        if cdf == 0 or Sm == 0:
                            cond_prob = pdf
                        else:
                            cond_prob = pdf/(get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1)) / Sm_dict[read_len]
                    elif config.LR_cond_prob_calc == 'form_2':
                        if cdf == 0:
                            cond_prob = pdf
                        else:
                            cond_prob = pdf/(get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1))
                    elif config.LR_cond_prob_calc == 'no_read_length':
                        cond_prob = 1
                    elif config.LR_cond_prob_calc == 'no_isoform_specific':
                        cond_prob = pdf
                if isoform_len - read_len + 1 != 0:
                    cond_prob /= isoform_len - read_len + 1
                else:
                    cond_prob = 0
                if isoform in isoform_index_dict:
                    isoform_index = isoform_index_dict[isoform]
                    if cond_prob == float('inf') or cond_prob <= MIN_PROB:
                        cond_prob = 0
                    all_reads_isoform_cond_prob[read_index,isoform_index] = cond_prob
        worker_cond_prob_matrix = convert_dict_to_sparse_matrix(all_reads_isoform_cond_prob,read_index+1,len(isoform_index_dict))
        all_cond_prob_matrix.append(worker_cond_prob_matrix)
        num_batches += 1
    cond_prob_matrix = scipy.sparse.vstack(all_cond_prob_matrix)
    Path(f'{output_path}/temp/cond_prob/').mkdir(exist_ok=True,parents=True)
    scipy.sparse.save_npz(f'{output_path}/temp/cond_prob/{worker_id}_cond_prob.npz',cond_prob_matrix)
    # print(f'Calculate the cond prob for LR: Worker {worker_id} done!',flush=True)
    return {worker_id:num_batches}
#     return all_reads_isoform_cond_prob
def get_cond_prob_MT_LIQA_modified(threads,output_path,isoform_df,isoform_index_dict,read_len_dist,Sm_dict):
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = worker_id,output_path,isoform_df,isoform_index_dict,read_len_dist,Sm_dict
        futures.append(pool.apply_async(get_all_reads_isoform_cond_prob_LIQA_modified,(args,)))
    num_batches_dict = {}
    for future in futures:
        num_batch = future.get()
        num_batches_dict.update(num_batch)
    pool.close()
    pool.join()
    return num_batches_dict

