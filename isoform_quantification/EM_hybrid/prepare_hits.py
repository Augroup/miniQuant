import pysam
from pathlib import Path
import pandas as pd
import pickle
import numpy as np
import multiprocessing as mp
import glob
import config
import shutil
import scipy
import scipy.sparse
import datetime
import time
from EM_hybrid.util import convert_dict_to_sparse_matrix,safe_divide_sparse
from EM_hybrid.cal_eff_len import get_eff_len_dict
def get_stats(arr):
    if len(arr) == 0:
        return np.array([float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')])
    s = scipy.stats.describe(arr)
    return np.array([s.minmax[0],s.minmax[1],s.mean,s.variance,s.skewness,s.kurtosis])
def dump_hits_dict(all_fragment_lengths,frag_len_dict,worker_id,batch_id,read_index,isoform_index_dict,output_path,theta_matrix):
# def dump_hits_dict(all_fragment_lengths,frag_len_dict,worker_id,batch_id,read_index,isoform_index_dict,output_path,theta_matrix,frag_len_sum,frag_len_squared_sum,num_frag_len):
    # st = time.time()
    chunk_frag_len_matrix = convert_dict_to_sparse_matrix(frag_len_dict,read_index,len(isoform_index_dict))
    # if frag_len_matrix is None:
    #     frag_len_matrix = chunk_frag_len_matrix
    # else:
    #     frag_len_matrix = scipy.sparse.vstack([frag_len_matrix,chunk_frag_len_matrix])
    compatib_matrix = chunk_frag_len_matrix.sign()
    # time1 = time.time() - st
    # st = time.time()
    # frag_len_matrix = convert_dict_to_sparse_matrix(frag_len_dict,read_index,len(isoform_index_dict))
    # compatib_dict = {}
    # for index in frag_len_dict:
    #     compatib_dict[index] = 1
    # compatib_matrix = convert_dict_to_sparse_matrix(compatib_dict,read_index,len(isoform_index_dict),dtype=np.int8)
    # compatib_matrix = csr_matrix(compatib_matrix)
    
    if config.inital_theta in ['SR_unique','LR_unique','hybrid_unique']:
        compatib_matrix = scipy.sparse.csr_matrix(compatib_matrix)
        unique_compatib_matrix = compatib_matrix[(compatib_matrix.sum(axis=1)==1).nonzero()[0]]
        chunk_theta_matrix = safe_divide_sparse(unique_compatib_matrix,unique_compatib_matrix.sum(axis=1)).sum(axis=0)
    elif config.inital_theta in ['SR','LR','hybrid','random','uniform']:
        compatib_matrix = scipy.sparse.csr_matrix(compatib_matrix)
        chunk_theta_matrix = safe_divide_sparse(compatib_matrix,compatib_matrix.sum(axis=1)).sum(axis=0)
    # time2 = time.time() - st
    # with open(f'{output_path}/temp/hits_dict/{worker_id}_hits_dict','ab') as f:
    #     pickle.dump([frag_len_dict,read_index,len(isoform_index_dict)],f)
    # theta_path = f'{output_path}/temp/hits_dict/{worker_id}_{batch_id}_chunk_theta_matrix.npz'
    # np.savez_compressed(theta_path,theta=chunk_theta_matrix)
    # st = time.time()
    if theta_matrix is None:
        theta_matrix = chunk_theta_matrix
    else:
        theta_matrix += chunk_theta_matrix
    # all_fragment_lengths_arr = np.array(all_fragment_lengths)
    # frag_len_sum += all_fragment_lengths_arr.sum()
    # frag_len_squared_sum += (all_fragment_lengths_arr*all_fragment_lengths_arr).sum()
    # num_frag_len += len(all_fragment_lengths)
    # with open(f'{output_path}/temp/fragment_lengths/{worker_id}_{batch_id}','wb') as f:
    #     pickle.dump(all_fragment_lengths,f)
    # frag_lengths_path = f'{output_path}/temp/fragment_lengths/{worker_id}_{batch_id}'
    # time3 = time.time() - st
    return theta_matrix,chunk_frag_len_matrix
    # ,all_fragment_lengths_arr,frag_len_sum,frag_len_squared_sum,num_frag_len,[time1,time2,time3]
def get_hits_dict(args):
    os.nice(10)
    worker_id,alignment_file_path,start_pos,end_pos,output_path,queue,isoform_index_dict,gene_isoform_index = args
    frag_len_dict = {}
    read_index = 0
    all_fragment_lengths = []
    previous_read_name = None
    theta_matrix = None
    all_frag_len_matrix = []
    all_frag_len_arr = []
    # times = [0,0,0]
    batch_id = 0
    buffer_size = 1e5
    # frag_len_sum,frag_len_squared_sum,num_frag_len = 0,0,0
    # with open(f'{output_path}/temp/hits_dict/{worker_id}_hits_dict','wb') as f:
    #     f.truncate(0)
    with pysam.AlignmentFile(alignment_file_path, "r",ignore_truncation=True) as f:
        # print(start_pos)
        f.seek(start_pos)
        for read in f:
            if read.template_length > 0:
                if '|' in read.reference_name:
                    isoform_name = read.reference_name.split('|')[0]
                else:
                    isoform_name = read.reference_name
                if isoform_name not in isoform_index_dict:
                    continue

                
                if previous_read_name is None:
                    previous_read_name = read.query_name
                elif previous_read_name != read.query_name:
                    if f.tell() > end_pos:
                        break
                    read_index += 1
                    previous_read_name = read.query_name
                    if len(frag_len_dict) >= buffer_size:
                        theta_matrix,frag_len_matrix = dump_hits_dict(all_fragment_lengths,frag_len_dict,worker_id,batch_id,read_index,isoform_index_dict,output_path,theta_matrix)
                        
                        # frag_len_arr,frag_len_sum,frag_len_squared_sum,num_frag_len,new_times = dump_hits_dict(all_fragment_lengths,frag_len_dict,worker_id,batch_id,read_index,isoform_index_dict,output_path,theta_matrix,frag_len_sum,frag_len_squared_sum,num_frag_len)
                        all_frag_len_matrix.append(frag_len_matrix)
                        all_frag_len_arr.append(np.array(all_fragment_lengths))
                        all_fragment_lengths = []
                        frag_len_dict = {}
                        batch_id += 1
                        read_index =0
            
                isoform_index = isoform_index_dict[isoform_name]
                frag_len_dict[(read_index,isoform_index)] = read.template_length
                # hits_dict[read.query_name].append({'fragment_length':read.template_length,'isoform':isoform_name})
                all_fragment_lengths.append(read.template_length)
    if len(frag_len_dict) > 0:
        read_index += 1
        theta_matrix,frag_len_matrix = dump_hits_dict(all_fragment_lengths,frag_len_dict,worker_id,batch_id,read_index,isoform_index_dict,output_path,theta_matrix)
        all_frag_len_matrix.append(frag_len_matrix)
        all_frag_len_arr.append(np.array(all_fragment_lengths))
        all_fragment_lengths = []
        frag_len_dict = {}
        batch_id += 1
        read_index =0
    if len(all_frag_len_arr) == 0:
        return {worker_id:batch_id},None,None
    # queue.put('done')
    frag_len_matrix = scipy.sparse.vstack(all_frag_len_matrix)
    theta_path = f'{output_path}/temp/hits_dict/{worker_id}_theta_matrix.npz'
    np.savez_compressed(theta_path,theta=theta_matrix)
    frag_len_path = f'{output_path}/temp/hits_dict/{worker_id}_frag_len_arr.npz'
    np.savez_compressed(frag_len_path,frag_len=np.hstack(all_frag_len_arr))
    frag_len_matrix_path = f'{output_path}/temp/hits_dict/{worker_id}_frag_len.npz'
    scipy.sparse.save_npz(frag_len_matrix_path,frag_len_matrix)

    hits_matrix = frag_len_matrix.sign()
    hits_matrix_sum = np.array(hits_matrix.sum(axis=1))[:,0]
    unique_mapping_mask = (hits_matrix_sum == 1)
    multi_mapping_mask = (hits_matrix_sum != 1)
    frag_len_matrix = scipy.sparse.csr_matrix(frag_len_matrix)
    unique_frag_len_matrix_csc = scipy.sparse.csc_matrix(frag_len_matrix[unique_mapping_mask,:])
    multi_frag_len_matrix_csc = scipy.sparse.csc_matrix(frag_len_matrix[multi_mapping_mask,:])
    SR_feature_dict = {}
    # duration_time0_list = []
    # duration_time1_list = []
    for gname,isoform_index in gene_isoform_index.items():
        # st = time.time()
        try:
            unique_mapping = unique_frag_len_matrix_csc[:,isoform_index].data
        except:
            unique_mapping = np.array([])
        # time0 = time.time()
        try:
            multi_mapping = multi_frag_len_matrix_csc[:,isoform_index].data
        except:
            multi_mapping = np.array([])
        # time1 = time.time()
        SR_feature_dict[gname] = {'unique_mapping':unique_mapping,'multi_mapping':multi_mapping}
        # duration_time0_list.append(time0-st)
        # duration_time1_list.append(time1-time0)
    # st = time.time()
    with open(f'{output_path}/temp/machine_learning/SR_feature_dict_{worker_id}','wb') as f:
        pickle.dump(SR_feature_dict,f)
    # dump_time = time.time() -st
    # print(f'Get_hits_dict: Worker {worker_id} done with {batch_id} batches!')
    # print(f'{worker_id}:'+str(len(duration_time0_list)))
    # print(f'{worker_id}:'+str(np.sum(duration_time0_list)))
    # print(f'{worker_id}:'+str(np.sum(duration_time1_list)))
    # print(f'{worker_id}:'+str(dump_time),flush=True)
    # print(f'Get_hits_dict: Worker {worker_id} done with {batch_id} batches!')
    # print(datetime.datetime.now(),flush=True)
    # print(times,flush=True)
    return {worker_id:batch_id},theta_path,frag_len_path
import os
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
    for i in range(threads):
        with open(alignment_file_path,'r') as f:
            if i == 0:
                start_offset = 0
                for line in f:
                    if line[0] != '@':
                        break
                    start_offset += len(line)
            else:
                f.seek(header_size+i*chunksize)
                f.readline()
                previous_read_name = None
                previous_byte_pos = f.tell()
                while True:
                    line = f.readline()
                    read_name = line.split('\t')[0]
                    if previous_read_name is None:
                        previous_read_name = read_name
                    elif previous_read_name != read_name:
                        start_offset = previous_byte_pos
                        break
                    previous_byte_pos = f.tell()
            byte_marker.append(start_offset)
    byte_marker.append(file_stats.st_size)
    return byte_marker
def get_all_hits_dict(alignment_file_path,byte_marker,threads,output_path,isoform_index_dict,gene_isoform_index):
    # watcher_pool = mp.Pool(1)
    worker_pool = mp.Pool(threads)
    manager = mp.Manager()
    queue = manager.Queue()    
    all_frag_len_arr_list = []
    # watcher = watcher_pool.apply_async(get_all_hits_dict_listener, args=(queue,threads,alignment_file_path,isoform_index_dict))
    futures = []
    num_batches_dict = {}
    for i in range(threads):
        start_pos,end_os = byte_marker[i],byte_marker[i+1]
        args = i,alignment_file_path,start_pos,end_os,output_path,queue,isoform_index_dict,gene_isoform_index
        futures.append(worker_pool.apply_async(get_hits_dict,(args,)))
    theta_arr = None
    for future in futures:
        # num_batch,theta_path,frag_len_path,worker_frag_len_sum,worker_frag_len_squared_sum,worker_num_frag_len = future.get()
        num_batch,theta_path,frag_len_path = future.get()
        if theta_path is None:
            continue
        # frag_len_sum += worker_frag_len_sum
        # frag_len_squared_sum += worker_frag_len_squared_sum
        # num_frag_len += worker_num_frag_len
        theta_matrix = np.load(theta_path)['theta']
        if theta_arr is None:
            theta_arr = theta_matrix
        else:
            theta_arr += theta_matrix
        worker_frag_len_arr = np.load(frag_len_path)['frag_len']
        all_frag_len_arr_list.append(worker_frag_len_arr)
        num_batches_dict.update(num_batch)
    # queue.put('kill')
    worker_pool.close()
    worker_pool.join()
    # watcher.get()
    all_frag_len_arr = np.hstack(all_frag_len_arr_list)
    mean_f_len = all_frag_len_arr.mean()
    std_f_len = all_frag_len_arr.std()
    # print(f'STD_F_LEN:{std_f_len}')
    # print(f'MEAN_F_LEN:{mean_f_len}')
    if config.mean_f_len is not None:
        mean_f_len = config.mean_f_len
    if config.std_f_len is not None:
        std_f_len = config.std_f_len
    eff_len_arr = get_eff_len_dict(alignment_file_path,mean_f_len,std_f_len,isoform_index_dict,config.eff_len_option)
    theta_arr = theta_arr/eff_len_arr
    config.pseudo_count_SR = 1/eff_len_arr.max()
    # theta_arr /= theta_arr.sum()
    # watcher_pool.close()
    # watcher_pool.join()
    return theta_arr,mean_f_len,std_f_len,eff_len_arr,num_batches_dict
def get_all_ant(args):
    worker_id,eff_len_arr,mean_f_len,std_f_len,output_path,num_batches = args
    num_batches_processed = 0
    frag_len_matrix_path = f'{output_path}/temp/hits_dict/{worker_id}_frag_len.npz'
    frag_len_matrix = scipy.sparse.load_npz(frag_len_matrix_path)
    ANT_matrix = scipy.sparse.csr_matrix(frag_len_matrix,dtype=float)
    ANT_matrix.data[:] = 1/(std_f_len * np.sqrt(np.pi)) * np.e**(-1/2*((ANT_matrix.data-mean_f_len)/std_f_len)**2)
    ANT_matrix = safe_divide_sparse(ANT_matrix,eff_len_arr)
    # try:
    #     Path(frag_len_matrix_path).unlink()
    # except:
    #     pass
    scipy.sparse.save_npz(f'{output_path}/temp/hits_dict/{worker_id}_ANT.npz',ANT_matrix)
    # print(f'Get_all_ant: Worker {worker_id} done!')
    # print(datetime.datetime.now(),flush=True)
def get_ant_all_workers(eff_len_arr,mean_f_len,std_f_len,threads,output_path,num_batches_dict,isoform_index_dict):
    pool = mp.Pool(threads)
    futures = []
    for i in range(threads):
        args = i,eff_len_arr,mean_f_len,std_f_len,output_path,num_batches_dict[i]
        futures.append(pool.apply_async(get_all_ant,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
def prepare_hits(SR_sam,output_path,isoform_index_dict,gene_isoform_index,threads):
    Path(f'{output_path}/temp/').mkdir(exist_ok=True,parents=True)
    # print('Sorting sam file by read name...',flush=True)
    # pysam.sort(SR_sam,'-n','-@',str(threads),'-o',f'{output_path}/temp/SR.sam')
    # alignment_file_path = f'{output_path}/temp/SR.sam'
    # pysam.collate(SR_sam,'-@',str(threads),'-o',f'{output_path}/temp/SR.sam')
    # alignment_file_path = f'{output_path}/temp/SR.sam'
    alignment_file_path = SR_sam
    # print('Done',flush=True)
    # print('Getting short reads info...',flush=True)
    # print(datetime.datetime.now(),flush=True)
    byte_marker = get_aln_line_marker(alignment_file_path,threads)
    Path(f'{output_path}/temp/hits_dict/').mkdir(exist_ok=True,parents=True)
    # Path(f'{output_path}/temp/fragment_lengths/').mkdir(exist_ok=True,parents=True)
    theta_arr,mean_f_len,std_f_len,eff_len_arr,num_batches_dict =  get_all_hits_dict(alignment_file_path,byte_marker,threads,output_path,isoform_index_dict,gene_isoform_index)
    # print('Get all hits dict done',flush=True)
    # print(datetime.datetime.now(),flush=True)
    get_ant_all_workers(eff_len_arr,mean_f_len,std_f_len,threads,output_path,num_batches_dict,isoform_index_dict)
    # print('Calculate ANT done',flush=True)
    # print(datetime.datetime.now(),flush=True)
    # try:
    #     Path(alignment_file_path).unlink()
    # except:
    #     pass
    return theta_arr,eff_len_arr,num_batches_dict
