from collections import defaultdict
import pysam
import time
import gc
import sys
import datetime
from pathlib import Path
from construct_feature_matrix import generate_all_feature_matrix_short_read
from construct_long_reads_feature_matrix import generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation, process_annotation_for_alignment
from parse_alignment_EM import parse_alignment_EM,get_region_read_count_length
from EM_libraries.EM_algo import EM_algo_main
from EM_theta_iterative_libraries.EM_algo import EM_algo_theta_iter_main
from EM_kde_libraries.EM_algo import EM_algo_kde_main
from EM_kde_score_libraries.EM_algo import EM_algo_kde_score_main
from EM_SR.EM_SR import EM_algo_SR
from EM_hybrid.EM_hybrid_community import EM_algo_hybrid
from EM_hybrid.EM_LR_alone_community import EM_algo_LR_alone
from construct_feature_matrix import calculate_all_condition_number,calculate_all_condition_number_long_read
import os
import config
from operator import itemgetter, attrgetter
import dill as pickle
import shutil
import pandas as pd
def infer_read_len(short_read_alignment_file_path):
    READ_LEN = 150
    sr_sam_valid = False
    with pysam.AlignmentFile(short_read_alignment_file_path, "r") as samfile:
            for read in samfile:
                READ_LEN = read.infer_query_length()
                if READ_LEN is not None:
                    sr_sam_valid = True
                    break
    if not sr_sam_valid:
        print('[ERROR] The SR sam seems invalid. Please check!')
        exit()
    return READ_LEN
def parse(ref_file_path, READ_JUNC_MIN_MAP_LEN, short_read_alignment_file_path, long_read_alignment_file_path,multi_mapping_filtering,threads):
    print('[INFO] Start calculating K-value...',flush=True)
    start_time = time.time()
    if short_read_alignment_file_path is not None:
        READ_LEN = infer_read_len(short_read_alignment_file_path)
    else:
        READ_LEN = 150
    gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict = parse_reference_annotation(
        ref_file_path, threads, READ_LEN, READ_JUNC_MIN_MAP_LEN, 'real_data')
    output_path =  config.output_path
    with open(f'{output_path}/temp/machine_learning/temp_LR_alignments_dict.pkl','wb') as f:
        pickle.dump([gene_isoforms_dict,LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict],f)
    short_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,150,allow_multi_exons=False) 
    SR_kvalue_dict = {}
    for rname in short_read_gene_matrix_dict:
        for gname in short_read_gene_matrix_dict[rname]:
            SR_kvalue_dict[gname] = short_read_gene_matrix_dict[rname][gname]['condition_number'][2]
    with open(f'{output_path}/temp/machine_learning/SR_kvalue_dict.pkl','wb') as f:
        pickle.dump(pd.Series(SR_kvalue_dict),f)
    _, gene_range, gene_interval_tree_dict = process_annotation_for_alignment(
        gene_exons_dict, gene_points_dict)
    end_time = time.time()
    print('[INFO] Done in %.3f s'%(end_time-start_time),flush=True)   
    return gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict, gene_range, gene_interval_tree_dict
# def map_short_reads(short_read_alignment_file_path, READ_JUNC_MIN_MAP_LEN, gene_isoforms_dict, gene_points_dict, gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, output_path, multi_mapping_filtering, threads):
#     print('Mapping short read to regions...', flush=True)
#     start_time = time.time()
#     if short_read_alignment_file_path is not None:
#         if multi_mapping_filtering == 'unique_only':
#             pysam.view('-F', '2816', '-q', '10', '-@', f'{threads}', '-h', '-o',
#                        f'{output_path}/temp_sr.sam', short_read_alignment_file_path, catch_stdout=False)
#             short_read_alignment_file_path = f'{output_path}/temp_sr.sam'
#         elif multi_mapping_filtering == 'best':
#             pysam.view('-F', '2816', '-@', f'{threads}', '-h', '-o',
#                        f'{output_path}/temp_sr.sam', short_read_alignment_file_path, catch_stdout=False)
#             short_read_alignment_file_path = f'{output_path}/temp_sr.sam'
#     short_read_gene_regions_read_count, SR_read_len, num_SRs = parse_alignment(short_read_alignment_file_path, READ_JUNC_MIN_MAP_LEN, gene_points_dict,
#                                                                                gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, False, False, threads)
#     config.READ_LEN = SR_read_len
#     print('Mapped {} short reads with read length {}'.format(
#         num_SRs, SR_read_len), flush=True)
#     end_time = time.time()
#     print('Constructing matrix and calculating condition number...', flush=True)
#     short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(
#         gene_isoforms_dict, SR_gene_regions_dict, short_read_gene_regions_read_count, SR_read_len, SR_genes_regions_len_dict, num_SRs)
#     print('Done in %.3f s' % (end_time-start_time), flush=True)
#     return short_read_gene_matrix_dict, SR_read_len
def map_long_reads(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,CHR_LIST,output_path,threads,multi_mapping_filtering):
    print('[INFO] Mapping long read to regions...',flush=True)
    start_time = time.time()
    if multi_mapping_filtering == 'unique_only':
        pysam.view('-F','2820','-q','10','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        pysam.sort(f'{output_path}/temp_lr.sam','-@',str(threads),'-m','3G','-o',f'{output_path}/temp_lr.sorted.sam')
        long_read_alignment_file_path = f'{output_path}/temp_lr.sorted.sam'
    elif multi_mapping_filtering == 'best':
        pysam.view('-F','2820','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        pysam.sort(f'{output_path}/temp_lr.sam','-@',str(threads),'-m','3G','-o',f'{output_path}/temp_lr.sorted.sam')
        long_read_alignment_file_path = f'{output_path}/temp_lr.sorted.sam'
    # print(long_read_alignment_file_path)
    parse_alignment_EM(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,output_path,threads,CHR_LIST)
    with open(f'{output_path}/temp/machine_learning/temp_LR_alignments_dict.pkl','rb') as f:
        [gene_isoforms_dict,LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict] = pickle.load(f)
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs = get_region_read_count_length(LR_gene_regions_dict,output_path) 
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,LR_gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,num_LRs,total_long_read_length,READ_JUNC_MIN_MAP_LEN,output_path,threads,False)
    LR_kvalue_dict = {}
    for rname in long_read_gene_matrix_dict:
        for gname in long_read_gene_matrix_dict[rname]:
            LR_kvalue_dict[gname] = long_read_gene_matrix_dict[rname][gname]['condition_number'][2]
    with open(f'{output_path}/temp/machine_learning/LR_kvalue_dict.pkl','wb') as f:
        pickle.dump(pd.Series(LR_kvalue_dict),f)
    # long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length,gene_regions_read_pos = parse_alignment(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
    # print('Mapped {} long reads'.format(num_LRs,flush=True))
    end_time = time.time()
    print('[INFO] Done in %.3f s'%(end_time-start_time),flush=True)   

import re
def parse_for_EM_algo(annotation):
    gene_exon_dict = {}
    gene_isoform_dict = {}
    isoform_exon_dict = {}
    isoform_gene_dict = {}
    strand_dict = {}
    with open(annotation,'r') as f:
        for line in f:
            if line.lstrip()[0] == "#":
                continue
            fields = line.split('\t')
            if (fields[2] != 'exon'):
                continue
            strand = fields[6]
            chr_name = fields[0]
            gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
            isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            if gene_name not in gene_exon_dict:
                gene_exon_dict[gene_name] = []
                gene_isoform_dict[gene_name] = set()
            if isoform_name not in isoform_exon_dict:
                isoform_exon_dict[isoform_name] = []
            gene_exon_dict[gene_name].append([start_pos,end_pos])
            gene_isoform_dict[gene_name].add(isoform_name)
            isoform_gene_dict[isoform_name] = gene_name
            isoform_exon_dict[isoform_name].append([start_pos,end_pos])
            if strand != '+' or strand != '-':
                strand = '+'
            strand_dict[gene_name] = strand
    for isoform in isoform_exon_dict:
        isoform_exon_dict[isoform] = sorted(isoform_exon_dict[isoform],key=lambda x:(x[0],x[1]))
    isoform_len_dict = {}
    for isoform in isoform_exon_dict:
        isoform_len_dict[isoform] = 0
        for exon in isoform_exon_dict[isoform]:
            isoform_len_dict[isoform] += exon[1] - exon[0] + 1
#     isoform_len_set = set(isoform_len_dict.values())
    return isoform_len_dict,isoform_exon_dict,strand_dict,isoform_gene_dict
def parse_and_dump_dict(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,multi_mapping_filtering,threads,READ_JUNC_MIN_MAP_LEN=15):
    _,gene_points_dict,gene_isoforms_dict,\
        _,_,LR_gene_regions_dict,LR_genes_regions_len_dict,\
            gene_isoforms_length_dict,raw_isoform_exons_dict,_,\
                _,_,gene_range,gene_interval_tree_dict = \
                    parse(ref_file_path,READ_JUNC_MIN_MAP_LEN,short_read_alignment_file_path,long_read_alignment_file_path,multi_mapping_filtering,threads)
    # check consistency of GTF and SAM
    gtf_isoforms = set()
    for rname in gene_isoforms_dict:
        for gname in gene_isoforms_dict[rname]:
            for isoform in gene_isoforms_dict[rname][gname]:
                gtf_isoforms.add(isoform)
    with pysam.AlignmentFile(short_read_alignment_file_path) as f:
        sam_isoforms = set(f.references)
    intersected = sam_isoforms.intersection(gtf_isoforms)
    intersected_length = len(intersected)
    # print(sam_isoforms)
    # print(gtf_isoforms)
    if intersected_length < len(sam_isoforms):
        if intersected_length == 0:
            print(f'[ERROR] There is no transcripts exist in both {short_read_alignment_file_path} and {ref_file_path}! Quiting now.')
            print(f"[ERROR] Please make sure to follow our recommended workflow in 'Data Preparation (if start from fasta/fastq file)' section of README.md to map the short reads")
            sys.exit(1)
        else:
            print(f'[WARNING] There are only {intersected_length} transcripts exist in both {short_read_alignment_file_path} and {ref_file_path}!')
            print(f"[WARNING] Please make sure to follow our recommended workflow in 'Data Preparation (if start from fasta/fastq file)' section of README.md to map the short reads")
    
    isoform_len_dict,isoform_exon_dict,strand_dict,isoform_gene_dict = parse_for_EM_algo(ref_file_path)

    start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST = dict(),dict(),dict(),dict(),list(gene_range.keys())
    CHR_LIST = list(gene_range.keys())
    for rname in CHR_LIST:     
        # Sort based on start position
        temp_list = sorted(gene_range[rname], key=itemgetter(1))
        start_pos_list[rname] = [temp_list[j][1] for j in range(len(temp_list))]
        start_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
        # Sort based on end position
        temp_list = sorted(gene_range[rname], key=itemgetter(2))
        end_pos_list[rname] = [temp_list[j][2] for j in range(len(temp_list))]
        end_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
    Path(f'{output_path}/temp/LR_alignments_dict/').mkdir(exist_ok=True,parents=True)
    for rname in gene_points_dict:
        dic = [gene_points_dict[rname],gene_interval_tree_dict[rname],LR_gene_regions_dict[rname],start_pos_list[rname],end_pos_list[rname],start_gname_list[rname],end_gname_list[rname]]
        with open(f'{output_path}/temp/LR_alignments_dict/{rname}','wb') as f:
            pickle.dump(dic,f)
    with open(f'{output_path}/temp/LR_alignments_dict/isoform_dict','wb') as f:
        pickle.dump([isoform_len_dict,isoform_exon_dict,strand_dict],f)
    with open(f'{output_path}/temp/LR_alignments_dict/gene_isoform_dict','wb') as f:
        pickle.dump([gene_isoforms_dict,isoform_len_dict,isoform_gene_dict],f)
    return CHR_LIST
            
    # for dic in [isoform_len_dict,isoform_exon_dict,strand_dict]:
        
    

def prepare_EM_LR(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,threads,multi_mapping_filtering='best',READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15,EM_choice='LIQA_modified',iter_theta='True'):
    CHR_LIST = parse_and_dump_dict(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,multi_mapping_filtering,threads)
    map_long_reads(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,CHR_LIST,output_path,threads,multi_mapping_filtering)
    with open(f'{output_path}/temp/LR_alignments_dict/gene_isoform_dict','rb') as f:
        [gene_isoforms_dict,isoform_len_dict,isoform_gene_dict] = pickle.load(f)
    try:
        shutil.rmtree(f'{output_path}/temp/LR_alignments_dict/')
    except:
        pass
    return isoform_len_dict,isoform_gene_dict,gene_isoforms_dict
def EM(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,alpha,beta,P,filtering,multi_mapping_filtering='best',SR_quantification_option='Mili',SR_fastq_list=[],reference_genome='',training=False,DL_model='',assign_unique_mapping_option='',threads=1,READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15,EM_choice='LIQA_modified',iter_theta='True'):
    # Path(output_path).mkdir(parents=True, exist_ok=True)
    # isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict = prepare_EM_LR(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,threads)
    # print(alpha)
    # print('Preprocessing...',flush=True)
    # print('Start quantification...',flush=True)
    # start_time = time.time()
    # if EM_choice == 'kde':
    #     EM_algo_kde_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice,config.kde_path,set(isoform_len_dict.values()))
    # elif EM_choice == 'kde_score':
    #     EM_algo_kde_score_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice,config.kde_path,set(isoform_len_dict.values()))
    # else:
    #     if iter_theta == 'True':
    #         EM_algo_theta_iter_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice)
    #     else:
    #         EM_algo_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice)
    # end_time = time.time()
    # print('Done in %.3f s'%(end_time-start_time),flush=True)
    pass
def EM_SR(short_read_alignment_file_path,output_path,threads):
    # Path(output_path).mkdir(parents=True, exist_ok=True)
    # start_time = time.time()
    # # isoform_len_dict,_,_ = parse_for_EM_algo(ref_file_path)
    # EM_algo_SR(short_read_alignment_file_path,output_path,threads)
    # end_time = time.time()
    # print('Done in %.3f s'%(end_time-start_time),flush=True)
    pass
def EM_hybrid(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,alpha,beta,P,filtering,multi_mapping_filtering='best',SR_quantification_option='Mili',SR_fastq_list=[],reference_genome='',training=False,DL_model='',assign_unique_mapping_option='',threads=1,READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15,EM_choice='LIQA_modified',iter_theta='True'):
    try:
        shutil.rmtree(output_path)
    except:
        pass
    Path(output_path).mkdir(parents=True, exist_ok=True)
    Path(f'{output_path}/temp/machine_learning/').mkdir(parents=True, exist_ok=True)
    file_stats = os.stat(long_read_alignment_file_path)
    total_bytes = file_stats.st_size
    # if total_bytes/1024/1024 < 100:
    #     threads = 1
    #     config.threads = 1
    print('[INFO] Preprocessing...',flush=True)
    isoform_len_dict,isoform_gene_dict,gene_isoforms_dict = prepare_EM_LR(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,threads)
    print('[INFO] Processing long reads data done at {}'.format(str(datetime.datetime.now())))
    print('[INFO] Start quantification...',flush=True)
    if short_read_alignment_file_path is None:
        EM_algo_LR_alone(isoform_len_dict,isoform_gene_dict,output_path,threads,EM_choice)
    else:
        file_stats = os.stat(short_read_alignment_file_path)
        total_bytes = file_stats.st_size
        # if total_bytes/1024/1024 < 100:
        #     threads = 1
        #     config.threads = 1
        EM_algo_hybrid(isoform_len_dict,isoform_gene_dict,gene_isoforms_dict,short_read_alignment_file_path,output_path,threads,EM_choice)
    try:
        Path(f'{output_path}/temp_sr.sam').unlink()
        Path(f'{output_path}/temp_sr.sorted.sam').unlink()
    except:
        pass
    try:
        Path(f'{output_path}/temp_lr.sam').unlink()
        Path(f'{output_path}/temp_lr.sorted.sam').unlink()
    except:
        pass
    try:
        shutil.rmtree(f'{output_path}/temp/')
    except:
        pass    
    try:
        Path(f'{output_path}/read_len_dist_sm_dict').unlink()
    except:
        pass    

