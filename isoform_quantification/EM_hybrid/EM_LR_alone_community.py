# hybrid
import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from pathlib import Path
import glob
import os
import gc
from EM_hybrid.EM_LR import prepare_LR
import scipy
import config
import datetime
import scipy.sparse
import time
import warnings
warnings.filterwarnings('ignore')

def E_step_LR(cond_prob,theta_arr):
    q = cond_prob.multiply(theta_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def M_step(isoform_q_arr_LR_all,theta_arr,eff_len_arr):
    ss = eff_len_arr / ((theta_arr * eff_len_arr).sum())
    # if alpha_df is None:
    new_theta_arr = ( isoform_q_arr_LR_all) / (isoform_q_arr_LR_all.sum())
    new_theta_arr = np.nan_to_num(new_theta_arr,0)
    # else:
    #     new_theta_arr = ((1-alpha_df) * isoform_q_arr_SR_all + alpha_df * isoform_q_arr_LR_all) / ((1-alpha_df) * isoform_q_df_SR.sum() * ss + alpha_df * isoform_q_arr_LR_all.sum())
    if new_theta_arr.sum() != 0:
        new_theta_arr = new_theta_arr/new_theta_arr.sum()
    new_theta_arr[new_theta_arr<1e-100] = 0
    return new_theta_arr
def build_dummy_gene_reads(isoform_gene_dict,isoform_index_dict,num_isoforms):
    gene_index_dict = {}
    for isoform,gene in isoform_gene_dict.items():
        if gene not in gene_index_dict:
            gene_index_dict[gene] = set()
        gene_index_dict[gene].add(isoform_index_dict[isoform])
    num_genes = 0
    row = []
    col = []
    data = []
    sorted_genes = sorted(gene_index_dict.keys())
    for gene in sorted_genes:
        isoform_indices = sorted(list(gene_index_dict[gene]))
        for isoform_index in isoform_indices:
            row.append(num_genes)
            col.append(isoform_index)
            data.append(1)
        num_genes += 1
    dummy_gene = scipy.sparse.coo_matrix((data, (row, col)), shape=(num_genes, num_isoforms))
    return dummy_gene
def get_connected_components(cond_prob_matrix):
    biadjacency = scipy.sparse.csr_matrix(cond_prob_matrix)
    adjacency = scipy.sparse.bmat([[None, biadjacency], [biadjacency.T, None]], format='csr')
    adjacency.sort_indices()
    n_components,labels = scipy.sparse.csgraph.connected_components(adjacency,return_labels=True)
    return labels
def serialize_community_worker(worker_id,worker_label,LR_labels,isoform_labels,cond_prob,output_path):
    worker_cond_prob = []
    worker_isoform = []
    worker_community = []
    for community_index in worker_label:
        community_LR_index = (LR_labels == community_index)
        community_isoforms_index = (isoform_labels == community_index)
        community_cond_prob = cond_prob[community_LR_index,:][:,community_isoforms_index]
        community_isoform = np.where(community_isoforms_index)[0]
        worker_cond_prob.append(community_cond_prob)
        worker_isoform.append(community_isoform)
        worker_community.append(community_index)
    np.savez(f'{output_path}/temp/community/{worker_id}.npz', cond_prob=worker_cond_prob,isoform=worker_isoform,community=worker_community)
def serialize_community(LR_labels,isoform_labels,cond_prob,output_path,threads):
    Path(f'{output_path}/temp/community/').mkdir(exist_ok=True,parents=True)
    pool = mp.Pool(threads)
    unique_labels = np.unique(isoform_labels)
    np.random.shuffle(unique_labels)
    all_worker_labels = np.array_split(unique_labels, threads)
    futures = []
    for worker_label,worker_id in zip(all_worker_labels,range(threads)):
       futures.append(pool.apply_async(serialize_community_worker,(worker_id,worker_label,LR_labels,isoform_labels,cond_prob,output_path,),error_callback=callback_error))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
def construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads):
    # ANT = scipy.sparse.vstack([scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{worker_id}_ANT.npz') for worker_id in range(threads)])
    cond_prob = scipy.sparse.vstack([scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{worker_id}_cond_prob.npz') for worker_id in range(threads)])
    # ANT = scipy.sparse.csr_matrix(ANT)
    cond_prob = scipy.sparse.csr_matrix(cond_prob)
    num_LRs,num_isoforms = cond_prob.shape[0],cond_prob.shape[1]
    # print(f'[INFO] Number of LRs:{num_LRs}')
    dummy_gene = build_dummy_gene_reads(isoform_gene_dict,isoform_index_dict,num_isoforms)
    cond_prob_matrix = scipy.sparse.vstack([cond_prob,dummy_gene])
    labels = get_connected_components(cond_prob_matrix)
    LR_labels,isoform_labels = labels[:num_LRs],labels[cond_prob_matrix.shape[0]:]
    serialize_community(LR_labels,isoform_labels,cond_prob,output_path,threads)
    return num_LRs
def EM_worker(worker_id,output_df,output_path,eff_len_arr,num_LRs):
    num_iters = config.EM_SR_num_iters
    min_diff = 1e-6
    worker_file = np.load(f'{output_path}/temp/community/{worker_id}.npz',allow_pickle=True)
    all_LR_expression_df = []
    all_community_iteration_df = []
    for community_cond_prob,community_isoform,community_id in zip(worker_file['cond_prob'],worker_file['isoform'],worker_file['community']):
        community_num_LRs = community_cond_prob.shape[0]
        theta_arr = np.ones(shape=(community_isoform.shape[0]))
        theta_arr = theta_arr/theta_arr.sum()
        community_eff_len_arr = np.array(eff_len_arr)[community_isoform]
        community_iteration_df = []
        for i in range(num_iters):
            theta_eff_len_product_arr = theta_arr * community_eff_len_arr /((theta_arr * community_eff_len_arr).sum())
            isoform_q_arr_LR = E_step_LR(community_cond_prob,theta_arr)
            new_theta_arr = M_step(isoform_q_arr_LR,theta_arr,community_eff_len_arr)
            new_theta_arr = np.array(new_theta_arr).flatten()
            diff = np.abs(theta_arr[new_theta_arr>1e-7] - new_theta_arr[new_theta_arr>1e-7])/new_theta_arr[new_theta_arr>1e-7]
            theta_arr = new_theta_arr
            # record iteration results
            iteration_df = pd.DataFrame({'theta':theta_arr},index=community_isoform)
            iteration_df['iteration'] = i
            community_iteration_df.append(iteration_df)
            if diff[diff > min_diff].shape[0] == 0:
                break
        # calculate LR TPM
        LR_expression = community_num_LRs/num_LRs * theta_arr * 1e6
        LR_expression_df = pd.DataFrame({'TPM':LR_expression,'theta':theta_arr},index=community_isoform)
        LR_expression_df['community'] = community_id
        LR_expression_df['community_num_LRs'] = community_num_LRs
        all_LR_expression_df.append(LR_expression_df)
        # write iteration results
        community_iteration_df = pd.concat(community_iteration_df)
        community_iteration_df = output_df.join(community_iteration_df,on='Index',how='inner').sort_values(['iteration','Isoform'])
        community_iteration_df['community'] = community_id
        all_community_iteration_df.append(community_iteration_df)
        # community_id = worker_file['community']
        # community_iteration_df.to_csv(f'{output_path}/EM_iterations/community_{community_id}.tsv',sep='\t',index=False)
    if len(all_LR_expression_df) == 0:
        LR_TPM_df = None
    else:
        all_LR_expression_df = pd.concat(all_LR_expression_df)
        LR_TPM_df = output_df.join(all_LR_expression_df,on='Index',how='inner')
    if len(all_community_iteration_df) == 0:
        all_community_iteration_df = None
    else: 
        all_community_iteration_df = pd.concat(all_community_iteration_df)
    return LR_TPM_df,all_community_iteration_df
def callback_error(result):
    print('ERR:', result,flush=True)
def EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads):
    print('[INFO] Start constructing the community...')
    st = time.time()
    num_LRs = construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads)
    duration = (time.time() - st)
    print(f'[INFO] Done in {duration} seconds!')
    print('[INFO] Start quantification...')  
    st = time.time()
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        futures.append(pool.apply_async(EM_worker,(worker_id,output_df,output_path,eff_len_arr,num_LRs,),error_callback=callback_error))
    all_LR_TPM_df = []
    all_iteration_df = []
    for future in futures:
        LR_TPM_df,iteration_df = future.get()
        all_LR_TPM_df.append(LR_TPM_df)
        all_iteration_df.append(iteration_df)
    # all_LR_TPM_df['TPM'] = all_LR_TPM_df['TPM']/all_LR_TPM_df['TPM'].sum() * 1e6
    pool.close()
    pool.join()
    all_LR_TPM_df = pd.concat(all_LR_TPM_df)
    all_LR_TPM_df['num_expected_LRs'] = all_LR_TPM_df['community_num_LRs'] * all_LR_TPM_df['theta']
    all_LR_TPM_df[['Isoform','Gene','TPM','theta','community','community_num_LRs','num_expected_LRs']].to_csv(f'{output_path}/LR_EM_expression.out',sep='\t',index=False)
    all_LR_TPM_df[['Isoform','Gene','TPM','num_expected_LRs']].sort_values(by=['Gene','Isoform']).to_csv(f'{output_path}/Isoform_abundance.out',sep='\t',index=False)
    all_iteration_df = pd.concat(all_iteration_df)
    all_iteration_df.to_csv(f'{output_path}/EM_iterations.tsv',sep='\t',index=False)
    duration = (time.time() - st)
    print(f'[INFO] Done in {duration} seconds!')
def EM_algo_LR_alone(isoform_len_dict,isoform_gene_dict,output_path,threads,EM_choice):
   # prepare arr
    isoform_len_df = pd.Series(isoform_len_dict)
    isoform_list = sorted(isoform_len_dict.keys())
    isoform_index_dict = {}
    isoform_len_arr = []
    for i,isoform in enumerate(isoform_list):
        isoform_index_dict[isoform] = i
        isoform_len_arr.append(isoform_len_dict[isoform])
    isoform_index_series = pd.Series(isoform_index_dict)
    isoform_index_series.name = 'Index'
    isoform_index_series.index.name = 'Isoform'
    output_df = isoform_index_series.to_frame().reset_index().sort_values(by="Index").set_index('Index')
    isoform_gene_df = pd.Series(isoform_gene_dict)
    isoform_gene_df.name = 'Gene'
    isoform_gene_df.index.name = 'Isoform'
    output_df = output_df.join(isoform_gene_df,on='Isoform')
    isoform_len_arr = np.array(isoform_len_arr)
    eff_len_arr = isoform_len_arr.copy()
    # prepare SR
    # theta_SR_arr,eff_len_arr,SR_num_batches_dict = prepare_hits(SR_sam,output_path,isoform_index_dict,threads)
    output_df['Effective length'] = eff_len_arr
    theta_LR_arr,_,LR_num_batches_dict = prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path)
    # num_SRs = theta_SR_arr.sum()
    num_LRs = theta_LR_arr.sum()
    # print(f'Number of SRs/eff_len:{num_SRs}')
    # print(f'Number of LRs:{num_LRs}')
    # print(f'Pseudo_count_SR:'+str(config.pseudo_count_SR))
    # print(f'Pseudo_count_LR:'+str(config.pseudo_count_LR),flush=True)
    # write_result_to_tsv(f'{output_path}/SR_count.tsv',output_df,theta_SR_arr.flatten())
    # write_result_to_tsv(f'{output_path}/LR_count.tsv',output_df,theta_LR_arr.flatten())
    # np.savez_compressed(f'{output_path}/initial_theta',theta=theta_arr)
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    # print('Using {} as initial theta'.format(config.inital_theta))
    
    # theta_arr = theta_arr/theta_arr.sum()
    EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads)
    


