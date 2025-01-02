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
from EM_hybrid.EM_SR import prepare_hits
from machine_learning.predict_alpha import predict_alpha
import scipy
import config
import datetime
import scipy.sparse
import scipy.stats
import time
import warnings
from EM_hybrid.util import sp_unique
warnings.filterwarnings('ignore')
def E_step_SR(ANT,theta_eff_len_product_arr):
    q = ANT.multiply(theta_eff_len_product_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def E_step_LR(cond_prob,theta_arr):
    q = cond_prob.multiply(theta_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def M_step(isoform_q_arr_SR_all,isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha):
    ss = eff_len_arr / ((theta_arr * eff_len_arr).sum())
    # if alpha_df is None:
    new_theta_arr = ((1-alpha) * isoform_q_arr_SR_all + alpha * isoform_q_arr_LR_all) / ((1-alpha) * isoform_q_arr_SR_all.sum() * ss + alpha * isoform_q_arr_LR_all.sum())
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
    gene_list = []
    for gene,isoform_indices in gene_index_dict.items():
        for isoform_index in isoform_indices:
            row.append(num_genes)
            col.append(isoform_index)
            data.append(1)
        num_genes += 1
        gene_list.append(gene)
    dummy_gene = scipy.sparse.coo_matrix((data, (row, col)), shape=(num_genes, num_isoforms))
    return dummy_gene,gene_list
def get_connected_components(cond_prob_matrix):
    biadjacency = scipy.sparse.csr_matrix(cond_prob_matrix)
    adjacency = scipy.sparse.bmat([[None, biadjacency], [biadjacency.T, None]], format='csr')
    adjacency.sort_indices()
    n_components,labels = scipy.sparse.csgraph.connected_components(adjacency,return_labels=True)
    return labels
def serialize_community_worker(worker_id,worker_label,isoform_labels,output_path):
    community_batch_index = worker_id
    worker_isoform = []
    worker_community = []
    for community_index in worker_label:
        community_isoforms_index = (isoform_labels == community_index)
        community_isoform = np.where(community_isoforms_index)[0]
        worker_isoform.append(community_isoform)
        worker_community.append(community_index)
    with open(f'{output_path}/temp/community/{community_batch_index}_isoform.pkl','wb') as f:
        pickle.dump(worker_isoform,f, protocol=4)
    with open(f'{output_path}/temp/community/{community_batch_index}_community.pkl','wb') as f:
        pickle.dump(worker_community,f, protocol=4)
    for reads_batch_id in range(config.threads):
        worker_ANT = []
        worker_cond_prob = []
        ANT = scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{reads_batch_id}_ANT.npz')
        ANT = scipy.sparse.csc_matrix(ANT)
        cond_prob = scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{reads_batch_id}_cond_prob.npz')
        cond_prob = scipy.sparse.csc_matrix(cond_prob)
        for community_index in worker_label:
            community_isoforms_index = (isoform_labels == community_index)
            community_ANT = scipy.sparse.csr_matrix(ANT[:,community_isoforms_index])
            # remove all zero rows
            community_ANT = community_ANT[community_ANT.getnnz(1)>0]
            community_cond_prob =  scipy.sparse.csr_matrix(cond_prob[:,community_isoforms_index])
            community_cond_prob = community_cond_prob[community_cond_prob.getnnz(1)>0]
            worker_ANT.append(community_ANT)
            worker_cond_prob.append(community_cond_prob)
        with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_ANT.pkl','wb') as f:
            pickle.dump(worker_ANT,f, protocol=4)
        with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_cond_prob.pkl','wb') as f:
            pickle.dump(worker_cond_prob,f, protocol=4)
    # for worker_label,community_batch_index in zip(all_worker_labels,range(len(all_worker_labels))):
    #     worker_ANT = []
    #     worker_cond_prob = []
    #     if worker_id == 0:
    #         worker_isoform = []
    #         worker_community = []
    #     for community_index in worker_label:
    #         community_SR_index = (SR_labels == community_index).nonzero()[0]
    #         community_LR_index = (LR_labels == community_index).nonzero()[0]
    #         community_isoforms_index = (isoform_labels == community_index)
    #         community_isoform = np.where(community_isoforms_index)[0]
    #         worker_community_SR_index = community_SR_index[(community_SR_index < num_processed_SR + ANT.shape[0]) & (community_SR_index >= num_processed_SR)].copy()
    #         worker_community_LR_index = community_LR_index[(community_LR_index < num_processed_LR + cond_prob.shape[0]) & (community_LR_index >= num_processed_LR)].copy()
    #         worker_community_SR_index -= num_processed_SR
    #         worker_community_LR_index -= num_processed_LR
    #         community_ANT = ANT[worker_community_SR_index,:][:,community_isoforms_index]
    #         community_cond_prob = cond_prob[worker_community_LR_index,:][:,community_isoforms_index]
    #         worker_ANT.append(community_ANT)
    #         worker_cond_prob.append(community_cond_prob)
    #         if worker_id == 0:
    #             worker_isoform.append(community_isoform)
    #             worker_community.append(community_index)
    #     with open(f'{output_path}/temp/community/{community_batch_index}_{worker_id}_ANT.pkl','wb') as f:
    #         pickle.dump(worker_ANT,f, protocol=4)
    #     with open(f'{output_path}/temp/community/{community_batch_index}_{worker_id}_cond_prob.pkl','wb') as f:
    #         pickle.dump(worker_cond_prob,f, protocol=4)
    #     if worker_id == 0:
    #         with open(f'{output_path}/temp/community/{community_batch_index}_isoform.pkl','wb') as f:
    #             pickle.dump(worker_isoform,f, protocol=4)
    #         with open(f'{output_path}/temp/community/{community_batch_index}_community.pkl','wb') as f:
    #             pickle.dump(worker_community,f, protocol=4)
def serialize_community(isoform_labels,output_path,threads):
    Path(f'{output_path}/temp/community/').mkdir(exist_ok=True,parents=True)
    pool = mp.Pool(threads)
    unique_labels = np.unique(isoform_labels)
    np.random.shuffle(unique_labels)
    all_worker_labels = np.array_split(unique_labels, threads)
    futures = []
    for worker_id in range(threads):
        futures.append(pool.apply_async(serialize_community_worker,(worker_id,all_worker_labels[worker_id],isoform_labels,output_path,),error_callback=callback_error))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
def construct_community_helper(isoform_gene_dict,isoform_index_dict,output_path,threads):
    st = time.time()
    ANT = []
    cond_prob = []
    # num_processed_SRs = []
    # num_processed_LRs = []
    for worker_id in range(threads):
        worker_ANT = scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{worker_id}_ANT.npz')
        worker_cond_prob = scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{worker_id}_cond_prob.npz')
        ANT.append(worker_ANT)
        cond_prob.append(worker_cond_prob)
        # num_processed_SRs.append(worker_ANT.shape[0])
        # num_processed_LRs.append(worker_cond_prob.shape[0])
    del worker_ANT
    del worker_cond_prob
    ANT = scipy.sparse.vstack(ANT)
    cond_prob = scipy.sparse.vstack(cond_prob)
    ANT = scipy.sparse.csr_matrix(ANT)
    cond_prob = scipy.sparse.csr_matrix(cond_prob)
    time1 = time.time()
    # if threads > 1:
    #     num_processed_SRs = [0]+[sum(np.array(num_processed_SRs)[:i]) for i in range(threads-1)]
    #     num_processed_LRs = [0]+[sum(np.array(num_processed_LRs)[:i]) for i in range(threads-1)]
    # else:
    #     num_processed_SRs = [0]
    #     num_processed_LRs = [0]
    num_SRs,num_LRs,num_isoforms = ANT.shape[0], cond_prob.shape[0],ANT.shape[1]
#     print(f'Number of SRs:{num_SRs}')
#     print(f'Number of LRs:{num_LRs}')
    dummy_gene,gene_list = build_dummy_gene_reads(isoform_gene_dict,isoform_index_dict,num_isoforms)
    time2 = time.time()
    reads_matrix_uniq = sp_unique(scipy.sparse.vstack([ANT.sign(),cond_prob.sign()]), axis=0)
    time3 = time.time()
    del ANT
    del cond_prob
    cond_prob_matrix = scipy.sparse.vstack([reads_matrix_uniq,dummy_gene])
    labels = get_connected_components(cond_prob_matrix)
    time4 = time.time()
    isoform_labels = labels[cond_prob_matrix.shape[0]:]
    gene_community_id_dict = {}
    for label,gene in zip(labels[reads_matrix_uniq.shape[0]:cond_prob_matrix.shape[0]],gene_list):
        gene_community_id_dict[gene] = label
    with open(f'{output_path}/temp/machine_learning/gene_community_id_dict.pkl','wb') as f:
        pickle.dump(gene_community_id_dict,f)
    # print(time1-st)
    # print(time2-time1)
    # print(time3-time2)
    # print(time4-time3)
    return isoform_labels,num_SRs,num_LRs
def construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads):
    isoform_labels,num_SRs,num_LRs = construct_community_helper(isoform_gene_dict,isoform_index_dict,output_path,threads)
    st = time.time()
    serialize_community(isoform_labels,output_path,threads)
    time5 = time.time()
    # print(time5-st)
    return num_SRs,num_LRs
def EM_worker(worker_id,output_df,output_path,eff_len_arr,num_SRs,num_LRs):
    num_iters = config.EM_SR_num_iters
    min_diff = 1e-6
    worker_ANT,worker_cond_prob,worker_isoform,worker_community = [],[],[],[]
    with open(f'{output_path}/temp/community/{worker_id}_isoform.pkl','rb') as f:
        worker_isoform = pickle.load(f)
    with open(f'{output_path}/temp/community/{worker_id}_community.pkl','rb') as f:
        worker_community = pickle.load(f)
    for i in range(config.threads):
        with open(f'{output_path}/temp/community/{worker_id}_{i}_ANT.pkl','rb') as f:
            temp_ANT = pickle.load(f)
        with open(f'{output_path}/temp/community/{worker_id}_{i}_cond_prob.pkl','rb') as f:
            temp_cond_prob = pickle.load(f)
        if i == 0:
            worker_ANT,worker_cond_prob = temp_ANT,temp_cond_prob
        else:
            for j in range(len(temp_ANT)):
                worker_ANT[j] = scipy.sparse.vstack([worker_ANT[j],temp_ANT[j]])
                worker_cond_prob[j] = scipy.sparse.vstack([worker_cond_prob[j],temp_cond_prob[j]])
    all_LR_expression_df = []
    all_SR_expression_df = []
    all_community_iteration_df = []
    for community_ANT,community_cond_prob,community_isoform,community_id in zip(worker_ANT,worker_cond_prob,worker_isoform,worker_community):
        if config.alpha_df_path is None:
            alpha = float(config.alpha)
        else:
            alpha_df = pd.read_csv(config.alpha_df_path,sep='\t')
            if len(alpha_df.columns) == 2:
                alpha_df['community_id'] = alpha_df['community_id'].astype(int)
                alpha_df = alpha_df.set_index('community_id')
                alpha = alpha_df.loc[community_id].values[0]
            elif len(alpha_df.columns) == 3:
                alpha_df['id'] = alpha_df['id'].astype(int)
                alpha_df = alpha_df.set_index('id')
                alpha = alpha_df.loc[community_id]['alpha']
        community_num_SRs = community_ANT.shape[0]
        community_num_LRs = community_cond_prob.shape[0]
        theta_arr = np.ones(shape=(community_isoform.shape[0]))
        theta_arr = theta_arr/theta_arr.sum()
        community_eff_len_arr = np.array(eff_len_arr)[community_isoform]
        community_iteration_df = []
        for i in range(num_iters):
            theta_eff_len_product_arr = theta_arr * community_eff_len_arr /((theta_arr * community_eff_len_arr).sum())
            isoform_q_arr_LR = E_step_LR(community_cond_prob,theta_arr)
            isoform_q_arr_SR = E_step_SR(community_ANT,theta_eff_len_product_arr)
            new_theta_arr = M_step(isoform_q_arr_SR,isoform_q_arr_LR,theta_arr,community_eff_len_arr,alpha)
            new_theta_arr = np.array(new_theta_arr).flatten()
            diff = np.abs(theta_arr[new_theta_arr>1e-7] - new_theta_arr[new_theta_arr>1e-7])/new_theta_arr[new_theta_arr>1e-7]
            # record iteration results
            iteration_df = pd.DataFrame({'theta':theta_arr},index=community_isoform)
            iteration_df['iteration'] = i
            community_iteration_df.append(iteration_df)
            theta_arr = new_theta_arr
            if diff[diff > min_diff].shape[0] == 0:
                break
        # calculate LR TPM
        LR_expression = community_num_LRs/num_LRs * theta_arr * 1e6
        LR_expression_df = pd.DataFrame({'TPM':LR_expression,'theta':theta_arr},index=community_isoform)
        LR_expression_df['community'] = community_id
        LR_expression_df['community_num_LRs'] = community_num_LRs
        all_LR_expression_df.append(LR_expression_df)
        # calculate SR TPM
        theta_eff_len_product_arr = theta_arr /((theta_arr * community_eff_len_arr).sum())
        theta_eff_len_product_arr = np.nan_to_num(theta_eff_len_product_arr,0)
        community_expression = (community_num_SRs * theta_eff_len_product_arr).sum()
        transcript_expression = community_num_SRs * theta_eff_len_product_arr
        SR_expression_df = pd.DataFrame({'transcript_expression':transcript_expression,'theta':theta_arr},index=community_isoform)
        SR_expression_df['community_expression'] = community_expression
        SR_expression_df['community'] = community_id
        SR_expression_df['community_num_SRs'] = community_num_SRs
        all_SR_expression_df.append(SR_expression_df)
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
    if len(all_SR_expression_df) == 0:
        SR_TPM_df = None
    else:
        all_SR_expression_df = pd.concat(all_SR_expression_df)
        SR_TPM_df = output_df.join(all_SR_expression_df,on='Index',how='inner')
    if len(all_community_iteration_df) == 0:
        all_community_iteration_df = None
    else: 
        all_community_iteration_df = pd.concat(all_community_iteration_df)
    return LR_TPM_df,SR_TPM_df,all_community_iteration_df
def callback_error(result):
    print('ERR:', result,flush=True)
def EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads,num_SRs,num_LRs):
    print('[INFO] Start quantification...')  
    st = time.time()
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        futures.append(pool.apply_async(EM_worker,(worker_id,output_df,output_path,eff_len_arr,num_SRs,num_LRs,),error_callback=callback_error))
    all_LR_TPM_df = []
    all_SR_TPM_df = []
    all_iteration_df = []
    for future in futures:
        LR_TPM_df,SR_TPM_df,iteration_df = future.get()
        all_LR_TPM_df.append(LR_TPM_df)
        all_SR_TPM_df.append(SR_TPM_df)
        all_iteration_df.append(iteration_df)
    # all_LR_TPM_df['TPM'] = all_LR_TPM_df['TPM']/all_LR_TPM_df['TPM'].sum() * 1e6
    pool.close()
    pool.join()
    all_LR_TPM_df = pd.concat(all_LR_TPM_df)
    all_LR_TPM_df['num_expected_LRs'] = all_LR_TPM_df['community_num_LRs'] * all_LR_TPM_df['theta']
    all_LR_TPM_df[['Isoform','Gene','TPM','theta','community','community_num_LRs']].to_csv(f'{output_path}/LR_EM_expression.out',sep='\t',index=False)
    all_SR_TPM_df = pd.concat(all_SR_TPM_df)
    all_SR_TPM_df['TPM'] = all_SR_TPM_df['community_expression']/(all_SR_TPM_df['transcript_expression'].sum()) * all_SR_TPM_df['theta'] * 1e6
    all_SR_TPM_df['num_expected_SRs'] = all_SR_TPM_df['community_num_SRs'] * all_SR_TPM_df['theta']
    all_SR_TPM_df[['Isoform','Gene','TPM','Effective length','theta','community','community_num_SRs']].to_csv(f'{output_path}/SR_EM_expression.out',sep='\t',index=False)
    all_SR_TPM_df = all_SR_TPM_df.set_index('Isoform').join(all_LR_TPM_df.set_index('Isoform')[['num_expected_LRs']]).reset_index()
    all_SR_TPM_df[['Isoform','Gene','Effective length','TPM','num_expected_SRs','num_expected_LRs']].sort_values(by=['Gene','Isoform']).to_csv(f'{output_path}/Isoform_abundance.out',sep='\t',index=False)
    all_iteration_df = pd.concat(all_iteration_df)
    all_iteration_df.to_csv(f'{output_path}/EM_iterations.tsv',sep='\t',index=False)
    duration = (time.time() - st)
    print('[INFO] Done in {} seconds at {}!'.format(duration,str(datetime.datetime.now())),flush=True)
def EM_algo_hybrid(isoform_len_dict,isoform_gene_dict,gene_isoforms_dict,SR_sam,output_path,threads,EM_choice):
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
    gene_isoform_index = {}
    for rname in gene_isoforms_dict:
        for gname in gene_isoforms_dict[rname]:
            gene_isoform_index[gname] = []
            for isoform in gene_isoforms_dict[rname][gname]:
                gene_isoform_index[gname].append(isoform_index_dict[isoform])
    print('[INFO] Start preparing short reads data... at {}'.format(str(datetime.datetime.now())),flush=True)
    theta_SR_arr,eff_len_arr,SR_num_batches_dict = prepare_hits(SR_sam,output_path,isoform_index_dict,gene_isoform_index,threads)
    output_df['Effective length'] = eff_len_arr
    print('[INFO] Prepare short reads data done at {}'.format(str(datetime.datetime.now())),flush=True)
    print('[INFO] Start preparing long reads data...',flush=True)
    theta_LR_arr,_,LR_num_batches_dict = prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path)
    num_SRs = theta_SR_arr.sum()
    num_LRs = theta_LR_arr.sum()
    print('[INFO] Prepare long reads data done at {}'.format(str(datetime.datetime.now())),flush=True)
    # print(f'Number of SRs/eff_len:{num_SRs}')
    # print(f'Number of LRs:{num_LRs}')
    # print(f'Pseudo_count_SR:'+str(config.pseudo_count_SR))
    # print(f'Pseudo_count_LR:'+str(config.pseudo_count_LR),flush=True)
    # write_result_to_tsv(f'{output_path}/SR_count.tsv',output_df,theta_SR_arr.flatten())
    # write_result_to_tsv(f'{output_path}/LR_count.tsv',output_df,theta_LR_arr.flatten())
    # np.savez_compressed(f'{output_path}/initial_theta',theta=theta_arr)
    # Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    # print('Using {} as initial theta'.format(config.inital_theta))
    
    # theta_arr = theta_arr/theta_arr.sum()
    print('[INFO] Start constructing the community...',flush=True)
    st = time.time()
    num_SRs,num_LRs = construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads)
    duration = (time.time() - st)
    print('[INFO] Done in {} seconds at {}!'.format(duration,str(datetime.datetime.now())),flush=True)
    print('[INFO] Extract features and predict best alpha...',flush=True)
    if config.alpha == 'adaptive':
        if not Path(config.alpha_df_path).exists():
            predict_alpha(output_path,num_SRs,num_LRs)
        else:
            print('[INFO] Using alpha from '+str(config.alpha_df_path))
    else:
        config.alpha_df_path = None
        print('[INFO] Using fixed alpha = '+str(config.alpha))
    EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads,num_SRs,num_LRs)
    


