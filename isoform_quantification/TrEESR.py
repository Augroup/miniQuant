from collections import defaultdict
from TransELS import infer_read_len
from construct_feature_matrix import calculate_all_condition_number,calculate_all_condition_number_long_read
from parse_annotation_main import parse_reference_annotation
from generate_output import generate_TrEESR_output
from construct_feature_matrix import generate_all_feature_matrix_short_read
from construct_long_reads_feature_matrix import generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment
from parse_alignment_main import parse_alignment
import config
import datetime
from pathlib import Path
import shutil
def TrEESR(ref_file_path,output_path,short_read_alignment_file_path,long_read_alignment_file_path,sr_region_selection,filtering,threads,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=0):
    try:
        shutil.rmtree(output_path)
    except:
        pass
    Path(output_path).mkdir(parents=True, exist_ok=True)
    Path(f'{output_path}/temp/machine_learning/').mkdir(parents=True, exist_ok=True)
    if short_read_alignment_file_path is not None:
        READ_LEN = infer_read_len(short_read_alignment_file_path)
    start_time = datetime.datetime.now()
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict,same_structure_isoform_dict,removed_gene_isoform_dict = \
        parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,sr_region_selection)
    end_time_1 = datetime.datetime.now()
    print('[INFO] Done in %.3f s'%((end_time_1-start_time).total_seconds()),flush=True)
    print('[INFO] Calculating the condition number...',flush=True)
    gene_regions_points_list,gene_range,gene_interval_tree_dict = process_annotation_for_alignment(gene_exons_dict,gene_points_dict)
    if short_read_alignment_file_path is not None:
        short_read_gene_regions_read_count, SR_read_len, num_SRs = parse_alignment(short_read_alignment_file_path, READ_JUNC_MIN_MAP_LEN, gene_points_dict,
                                                                               gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, False, False, threads)
        config.READ_LEN = SR_read_len
        short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(gene_isoforms_dict, SR_gene_regions_dict, short_read_gene_regions_read_count, SR_read_len, SR_genes_regions_len_dict, num_SRs,False)
    else:
        SR_read_len = 150
        short_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,SR_read_len,allow_multi_exons=False)
    if long_read_alignment_file_path is not None:
        long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length,_ = parse_alignment(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
        long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,LR_gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,num_LRs,total_long_read_length,READ_JUNC_MIN_MAP_LEN,output_path,threads,False)
    else:
        long_read_gene_matrix_dict = calculate_all_condition_number_long_read(gene_isoforms_dict,LR_gene_regions_dict,allow_multi_exons=True)
    end_time_2 = datetime.datetime.now()
    print('[INFO] Done in %.3f s'%((end_time_2-end_time_1).total_seconds()),flush=True)
    raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict = defaultdict(dict),defaultdict(dict),defaultdict(dict)
    raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            raw_gene_num_exon_dict[chr_name][gene_name] = len(raw_gene_exons_dict[chr_name][gene_name])
            gene_num_exon_dict[chr_name][gene_name] = len(gene_exons_dict[chr_name][gene_name])
            gene_num_isoform_dict[chr_name][gene_name] = len(gene_isoforms_dict[chr_name][gene_name])
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name]:
                raw_isoform_num_exon_dict[isoform_name] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                isoform_length_dict[isoform_name] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                num_isoforms_dict[isoform_name] =  len(raw_isoform_exons_dict[chr_name][gene_name])
    info_dict_list = [raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict,raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict]
    gene_feature_dict = generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list,same_structure_isoform_dict,removed_gene_isoform_dict,gene_points_dict)
    ref_annotation_dict_list = [gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict]
    return gene_feature_dict
def get_kvalues_dict(ref_file_path,threads,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=10):
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=True)
    return long_read_gene_matrix_dict


