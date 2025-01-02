import argparse
from TrEESR import TrEESR
from TransELS import TransELS
from EM import EM,EM_SR,EM_hybrid
import config
import os
# import os
# os.system("taskset -p 0xfffff %d" % os.getpid())
# affinity_mask = os.sched_getaffinity(0)
# os.sched_setaffinity(0, affinity_mask)
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_TrEESR = subparsers.add_parser('cal_K_value', aliases=['TrEESR'],help='Calculate K values')
    # parser_TransELS = subparsers.add_parser('quantify', aliases=['TransELS'],help='Isoform quantification')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for calculation of K value')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=False)
    optional_TrEESR.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False)
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_TrEESR.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:read_length][read_length,num_exons,real_data]")
    optional_TrEESR.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")
    optional_TrEESR.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:True][True,False]")
    optional_TrEESR.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_TrEESR.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_TrEESR.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_TrEESR.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_TrEESR.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_TrEESR.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_TrEESR.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_TrEESR.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    weight_path = os.path.dirname(os.path.realpath(__file__))+'/weights/nanosim_weight_dict.pkl'
    # assert os.path.exists(weight_path)
    optional_TrEESR.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")

    # requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for isoform quantification')
    # requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    # requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    # requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    # optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    # optional_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    # optional_TransELS.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    # optional_TransELS.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    # optional_TransELS.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    # optional_TransELS.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    # optional_TransELS.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, kallisto,Salmon, RSEM] [default:kallisto]",default='kallisto')
    # # optional_TransELS.add_argument('--kallisto_index', type=str, help="kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/kallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    # optional_TransELS.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    # optional_TransELS.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    # optional_TransELS.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    # optional_TransELS.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    # optional_TransELS.add_argument('--training',type=str,default='False', help="Generate training dict")
    # optional_TransELS.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    # optional_TransELS.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    # optional_TransELS.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    # optional_TransELS.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    # optional_TransELS.add_argument('--use_weight_matrix',type=str, default='True',help="Whether use weight matrix[default:True][True,False]")
    # optional_TransELS.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    # optional_TransELS.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    # optional_TransELS.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    # optional_TransELS.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    # optional_TransELS.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    # optional_TransELS.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    # optional_TransELS.add_argument('--sr_region_selection',type=str, default='real_data',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    # optional_TransELS.add_argument('--keep_sr_exon_region',type=str, default='True',help="Keep exon region for SR if using real data to filter region [default:True][True,False]")
    # optional_TransELS.add_argument('--region_weight_path',type=str, default=weight_path,help="Mili LR region weight path")
   
    parser_EM = subparsers.add_parser('quantify', aliases=['EM'],help='Isoform quantification by EM algorithm')
    requiredNamed_EM = parser_EM.add_argument_group('required named arguments for isoform quantification')
    requiredNamed_EM.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_EM.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_EM = parser_EM.add_argument_group('optional arguments')
    optional_EM.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False,default=None)
    optional_EM.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    optional_EM.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    optional_EM.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    optional_EM.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    optional_EM.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    optional_EM.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, kallisto,Salmon, RSEM] [default:kallisto]",default='kallisto')
    # optional_EM.add_argument('--kallisto_index', type=str, help="kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/kallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    optional_EM.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    optional_EM.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    optional_EM.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    optional_EM.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    optional_EM.add_argument('--training',type=str,default='False', help="Generate training dict")
    optional_EM.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    optional_EM.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    optional_EM.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_EM.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_EM.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_EM.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    # optional_EM.add_argument('--same_struc_isoform_handling',type=str, default='keep',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_EM.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_EM.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_EM.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    optional_EM.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_EM.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    optional_EM.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    optional_EM.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_EM.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")
    optional_EM.add_argument('--EM_choice',type=str, default='LR',help="EM_choice[SR,LR,hybrid]")
    optional_EM.add_argument('--iter_theta',type=str, default='False',help="Whether use updated theta to re-calculate conditional prob [True,False]")
    optional_EM.add_argument('--kde_path',type=str, default='/fs/project/PCON0009/Au-scratch2/haoran/_projects/long_reads_rna_seq_simulator/models/kde_H1-hESC_dRNA',help="KDE model path")
    optional_EM.add_argument('--eff_len_option',type=str, default='kallisto',help="Calculation of effective length option [kallisto,RSEM]")
    optional_EM.add_argument('--EM_SR_num_iters',type=int, default=200,help="Number of EM SR iterations")
    optional_EM.add_argument('--EM_output_frequency',type=int, default=200,help="Frequency(in itertations) of outputting EM results")
    optional_EM.add_argument('--pretrained_model_path',type=str, default='cDNA-ONT',help="The pretrained model path to identify the alpha")
    optional_EM.add_argument('--alpha_df_path',type=str, default=None,help="Alpha df path")
    optional_EM.add_argument('--inital_theta','--initial_theta',type=str, default='uniform',help="initial_theta [LR,SR,LR_unique,SR_unique,uniform,hybrid,hybrid_unique,random]")
    optional_EM.add_argument('--inital_theta_eps','--initial_theta_eps',type=float, default=0.0,help="initial_theta eps [float]")
    optional_EM.add_argument('--eps_strategy',type=str, default='add_eps_small',help="how to add initial_theta eps [add_eps_all,add_eps_small]. (add_eps_small: add isoform with theta < eps with eps. add_eps: add eps to all isoforms)")
    optional_EM.add_argument('--isoform_start_end_site_tolerance',type=int, default=20,help="Isoform Start and end site tolerance for mapping long reads")
    optional_EM.add_argument('--junction_site_tolerance',type=int, default=5,help="Junction site tolerance for mapping long reads")
    optional_EM.add_argument('--read_len_dist_sm_dict_path',type=str, default=None,help="The path of read length distribution for long reads")
    optional_EM.add_argument('--LR_cond_prob_calc',type=str, default='form_2',help="How to calculate LR length distribution [form_1,form_2]")
    optional_EM.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")

    args = parser.parse_args()
    if args.filtering == 'True':
        args.filtering = True
    else:
        args.filtering = False
    # config.same_struc_isoform_handling = args.same_struc_isoform_handling
    config.output_path = args.output_path
    config.threads = args.threads
    config.same_struc_isoform_handling = 'keep'
    config.READ_JUNC_MIN_MAP_LEN = args.READ_JUNC_MIN_MAP_LEN
    config.multi_exon_region_weight = args.multi_exon_region_weight
    config.sr_region_selection = args.sr_region_selection
    config.region_weight_path = args.region_weight_path
    config.sr_design_matrix = args.sr_design_matrix
    if args.output_matrix_info == 'True':
        config.output_matrix_info = True
    else:
        config.output_matrix_info = False
    config.keep_sr_exon_region = args.keep_sr_exon_region
    if args.normalize_sr_A == 'True':
        config.normalize_sr_A = True
    else:
        config.normalize_sr_A = False
    if args.normalize_lr_A == 'True':
        config.normalize_lr_A = True
    else:
        config.normalize_lr_A = False
    if args.use_weight_matrix == 'True':
        config.use_weight_matrix = True
    else:
        config.use_weight_matrix = False
    config.add_full_length_region = args.add_full_length_region
    config.singular_values_tol = args.singular_values_tol
    # config.kallisto_index = args.kallisto_index
    # print('\n'.join(f'{k}={v}' for k, v in vars(args).items()))
    if args.subparser_name in ['cal_K_value','TrEESR']:
        print('[INFO] Calculate K values')
        TrEESR(args.gtf_annotation_path,args.output_path,args.short_read_sam_path,args.long_read_sam_path,args.sr_region_selection,args.filtering,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    # elif args.subparser_name in ['quantify','TransELS']:
    #     if args.training == 'True':
    #         args.training = True
    #     else:
    #         args.training = False
    #     print('Isoform quantification',flush=True)
    #     if (args.short_read_sam_path is None) or (args.alpha == 1.0):
    #         args.alpha = 1.0
    #         args.SR_quantification_option = 'Mili'
    #     if args.alpha != 1.0:
    #         if args.short_read_sam_path is None:
    #             raise Exception('You need to provide a short read alignemnt file if the alpha is not 1!')
    #         if args.SR_quantification_option != 'Mili':
    #             if (args.short_read_fastq is None) and (args.short_read_mate1_fastq is None or args.short_read_mate2_fastq is None):
    #                 raise Exception('You need to provide the single end or paried end SR fastq if using other short read quantification options!')
    #             if (args.reference_genome is None):
    #                 raise Exception('You need to provide the reference genome if using other short read quantification options!')
    #     if (args.alpha == 'adaptive'):
    #         alpha = 'adaptive'
    #     else:
    #         try:
    #             alpha = float(args.alpha)
    #         except:
    #             raise Exception('Alpha given is not numeric')
    #     if (args.beta == 'adaptive'):
    #         beta = 'adaptive'
    #     else:
    #         try:
    #             beta = float(args.beta)
    #         except:
    #             raise Exception('Beta given is not numeric')
    #     if args.SR_quantification_option not in ['Mili','kallisto','Salmon','RSEM']:
    #         raise Exception('SR_quantification_option is not valid.Options: [Mili, kallisto,Salmon, RSEM]')
    #     if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
    #         args.multi_mapping_filtering = 'no_filtering'
    #     SR_fastq_list = []
    #     if args.short_read_fastq is not None:
    #         SR_fastq_list = [args.short_read_fastq]
    #     elif args.short_read_mate1_fastq is not None:
    #         SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
    #     if args.DL_model is None:
    #         args.DL_model = args.SR_quantification_option + '.pt'
    #     TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    elif args.subparser_name in ['quantify','EM']:
        config.kde_path = args.kde_path
        if args.training == 'True':
            args.training = True
        else:
            args.training = False
        print('[INFO] Isoform quantification by miniQuant',flush=True)
        if (args.short_read_sam_path is None) or (args.alpha == 1.0):
            args.alpha = 1.0
            args.SR_quantification_option = 'Mini'
        # if args.alpha != 1.0:
        #     if args.short_read_sam_path is None:
        #         raise Exception('You need to provide a short read alignemnt file if the alpha is not 1!')
        #     if args.SR_quantification_option != 'Mili':
        #         if (args.short_read_fastq is None) and (args.short_read_mate1_fastq is None or args.short_read_mate2_fastq is None):
        #             raise Exception('You need to provide the single end or paried end SR fastq if using other short read quantification options!')
        #         if (args.reference_genome is None):
        #             raise Exception('You need to provide the reference genome if using other short read quantification options!')
        if (args.alpha == 'adaptive'):
            alpha = 'adaptive'
        else:
            try:
                alpha = float(args.alpha)
            except:
                raise Exception('Alpha given is not numeric')
        if (args.beta == 'adaptive'):
            beta = 'adaptive'
        else:
            try:
                beta = float(args.beta)
            except:
                raise Exception('Beta given is not numeric')
        # if args.SR_quantification_option not in ['Mili','kallisto','Salmon','RSEM']:
        #     raise Exception('SR_quantification_option is not valid.Options: [Mili, kallisto,Salmon, RSEM]')
        if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
            args.multi_mapping_filtering = 'no_filtering'
        SR_fastq_list = []
        if args.short_read_fastq is not None:
            SR_fastq_list = [args.short_read_fastq]
        elif args.short_read_mate1_fastq is not None:
            SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
        if args.DL_model is None:
            args.DL_model = args.SR_quantification_option + '.pt'
        config.EM_SR_num_iters = args.EM_SR_num_iters
        config.inital_theta_eps = args.inital_theta_eps
        config.EM_output_frequency = args.EM_output_frequency
        config.isoform_start_end_site_tolerance = args.isoform_start_end_site_tolerance
        config.junction_site_tolerance = args.junction_site_tolerance
        config.eps_strategy = args.eps_strategy
        config.read_len_dist_sm_dict_path = args.read_len_dist_sm_dict_path
        config.LR_cond_prob_calc = args.LR_cond_prob_calc
        if args.pretrained_model_path in ['cDNA-ONT','dRNA-ONT','cDNA-PacBio']:
            args.pretrained_model_path = os.path.dirname(os.path.realpath(__file__))+'/pretrained_models/' + args.pretrained_model_path +'/'
        config.pretrained_model_path = args.pretrained_model_path
        if args.EM_choice == 'SR':
            config.eff_len_option = args.eff_len_option
            args.long_read_sam_path = None
            args.alpha = 0
            args.inital_theta = 'SR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        elif args.EM_choice == 'hybrid':
            # args.alpha = 0.5
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            if args.alpha_df_path is None:
                config.alpha_df_path = args.output_path +'/hybrid_alpha.tsv'
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        else:
            if args.EM_choice == 'LR':
                args.EM_choice = 'LIQA_modified'
            args.short_read_sam_path = None
            args.alpha = 1
            args.inital_theta = 'LR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)    
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()
