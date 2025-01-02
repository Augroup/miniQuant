from collections import defaultdict
import numpy as np
from numpy import linalg as LA
# from qpsolvers import solve_qp
# from predict_params import predict_params_all_genes
import config
# from tqdm import tqdm
def normalize_expression(gene_isoform_expression_dict):
    gene_isoform_tpm_expression_dict = defaultdict(lambda: defaultdict(dict))
    # SR_isoform_expression_sum = 0
    LR_isoform_expression_sum = 0
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            # gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_expected_counts'] = gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts']
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_expected_counts'] = gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression']
            # SR_isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'].sum()
            LR_isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'].sum()

    # if SR_isoform_expression_sum == 0:
    #     SR_isoform_expression_sum = 1
    if LR_isoform_expression_sum == 0:
        LR_isoform_expression_sum = 1 
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_tpm'] = gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression']
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_tpm'] = gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'] * 1e6 / LR_isoform_expression_sum
            gene_alpha = gene_isoform_expression_dict[chr_name][gene_name]['alpha']
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['alpha'] = gene_alpha
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'] = (1 - gene_alpha) * gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] + gene_alpha * gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_tpm']

    return gene_isoform_tpm_expression_dict 
                

def estimate_isoform_expression_grid_search_iteration(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P,assign_unique_mapping_option,SR_quantification_option,params):
    if SR_region_read_count_matrix.sum() != 0:
        SR_b = SR_region_read_count_matrix / SR_region_read_count_matrix.sum()
    else:
        SR_b = SR_region_read_count_matrix.copy()
    # if LR_region_read_count_matrix.sum() != 0:
    #     LR_b = LR_region_read_count_matrix / LR_region_read_count_matrix.sum()
    # else:
    #     LR_b = LR_region_read_count_matrix.copy()
    LR_b = LR_region_read_count_matrix.copy()
    alpha,beta = params['alpha'], params['beta']
    if SR_quantification_option != 'Mili':
        alpha = 1.0
    num_isoforms = SR_isoform_region_matrix.shape[1]
    #l2 norm
    Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta * np.identity(num_isoforms)
    c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_b.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_b.T)
    lb = np.zeros(num_isoforms)
    if alpha == 1.0:
        G = - LR_isoform_region_matrix[np.count_nonzero(LR_isoform_region_matrix,axis=1) == 1,:]
        h = - LR_b[np.count_nonzero(LR_isoform_region_matrix,axis=1) == 1]
    elif alpha == 0.0:
        G = - SR_isoform_region_matrix[np.count_nonzero(SR_isoform_region_matrix,axis=1) == 1,:]
        h = - SR_b[np.count_nonzero(SR_isoform_region_matrix,axis=1) == 1]
    else:
        G = - np.concatenate((SR_isoform_region_matrix[np.count_nonzero(SR_isoform_region_matrix,axis=1) == 1,:], LR_isoform_region_matrix[np.count_nonzero(LR_isoform_region_matrix,axis=1) == 1,:]), axis=0)
        h = - np.concatenate((SR_b[np.count_nonzero(SR_isoform_region_matrix,axis=1) == 1], LR_b[np.count_nonzero(LR_isoform_region_matrix,axis=1) == 1]), axis=0)
    # l1 norm
    # Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix)
    # c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T) + beta* np.ones(num_isoforms)
    # lb = np.zeros(num_isoforms)
    # G = - np.concatenate((SR_isoform_region_matrix[SR_region_read_count_matrix>0,:], LR_isoform_region_matrix[LR_region_read_count_matrix>0,:]), axis=0)
    # h = - np.ones(G.shape[0])/(1/P)
    if assign_unique_mapping_option == 'linear_model':
        if config.use_weight_matrix:
            isoform_expression = solve_qp(Q, c,lb = lb)
        else:
            isoform_expression = solve_qp(Q, c,G,h, lb = lb)
    else:
        isoform_expression = solve_qp(Q, c,lb = lb)
    # isoform_expression = solve_qp(Q, c,lb = lb)
    # if ((isoform_expression+1e-6<0).any()):
    #     raise ValueError('Obtain negative value for isoform expression')
    isoform_expression[isoform_expression<0] = 0
    # target = (1.0 - alpha) * LA.norm(SR_region_read_count_matrix - np.matmul(SR_isoform_region_matrix,isoform_expression)) + alpha * LA.norm(LR_region_read_count_matrix - np.matmul(LR_isoform_region_matrix,isoform_expression))
    # 
    # perfect_isoform_expression = np.matmul(LA.inv((1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix)+ beta * np.identity(num_isoforms)), (1 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix) )
    # return isoform_expression,perfect_isoform_expression,target
    return isoform_expression

def assign_reads(isoform_expression,A,b):
    # unique mapping reads
    # A = A.copy()
    unique_regions_index = np.count_nonzero(A,axis=1) == 1
    unique_A = A[unique_regions_index]
    unique_b = b[unique_regions_index]
    unique_A[unique_A!=0] = 1
    unique_isoform_expression = unique_b.dot(unique_A)
    # multi mapping reads
    multi_regions_index = (np.count_nonzero(A,axis=1) > 1) & (b!=0)
    multi_A = A[multi_regions_index]
    multi_b = b[multi_regions_index]
    if (multi_A.dot(isoform_expression) == 0).any():
        if isoform_expression.sum() == 0:
            theta = isoform_expression + 1
            theta = theta/theta.sum()
        else:
            theta = isoform_expression + 0.01 * isoform_expression[isoform_expression>0].min()
    else:
        theta = isoform_expression
    X = multi_b/(multi_A.dot(theta))
    Y = theta * multi_A
    multiple_isoform_expression = X.dot(Y)
    corrected_isoform_expression = unique_isoform_expression + multiple_isoform_expression
    return corrected_isoform_expression
def estimate_isoform_expression_single_gene(args):
    
    short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,alpha,beta,P,model,assign_unique_mapping_option,SR_quantification_option,gene = args
    SR_isoform_region_matrix = short_read_gene_matrix_dict['isoform_region_matrix']
    SR_region_read_count_matrix = short_read_gene_matrix_dict['region_abund_matrix']
    SR_gene_counts = short_read_gene_matrix_dict['num_SRs_mapped_gene']
    SR_region_eff_length_matrix = np.zeros(SR_region_read_count_matrix.shape[0])
    for region_name  in short_read_gene_matrix_dict['region_names_indics']:
        index = short_read_gene_matrix_dict['region_names_indics'][region_name]
        SR_region_eff_length_matrix[index] = short_read_gene_matrix_dict['region_eff_length_dict'][region_name]
        assert short_read_gene_matrix_dict['region_eff_length_dict'][region_name] != 0
    LR_isoform_region_matrix = long_read_gene_matrix_dict['isoform_region_matrix']
    LR_region_read_count_matrix = long_read_gene_matrix_dict['region_abund_matrix']
    isoform_lengths = np.zeros((len(gene_isoforms_length_dict)))
    isoform_num_exons = np.zeros((len(gene_isoforms_length_dict)))
    isoform_names_indics = short_read_gene_matrix_dict['isoform_names_indics']
    for isoform_name in isoform_names_indics:
        isoform_lengths[isoform_names_indics[isoform_name]] = gene_isoforms_length_dict[isoform_name]
        isoform_num_exons[isoform_names_indics[isoform_name]] = long_read_gene_matrix_dict['num_exons'][isoform_name]
    num_isoforms = SR_isoform_region_matrix.shape[1]
    if ((SR_region_read_count_matrix<=0).all() and (LR_region_read_count_matrix<=0).all()):
        return np.zeros(num_isoforms),1.0
    if (beta == 'adaptive'):
        # beta_selections = [10**(-i) for i in range(1,10)]
        selected_beta = 1e-6
    else:
        selected_beta = beta
    # use alpha = 1.0 for Mili + other SR methods for linear model
    params = {'alpha':1.0,'beta':selected_beta}
    isoform_expression = estimate_isoform_expression_grid_search_iteration(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P,assign_unique_mapping_option,SR_quantification_option,params)
    if isoform_expression.sum() != 0:
        isoform_expression = isoform_expression/isoform_expression.sum()
    if assign_unique_mapping_option == 'linear_model':
        # SR_expression = isoform_lengths * isoform_expression
        # if SR_expression.sum() != 0:
        #     SR_expression = SR_expression/SR_expression.sum()
        # SR_expected_counts = SR_gene_counts * SR_expression
        # SR_isoform_expression = SR_expected_counts / isoform_lengths
        LR_isoform_expression = LR_region_read_count_matrix.sum() * isoform_expression
    else:
        # SR_expected_counts = assign_reads(isoform_expression,SR_isoform_region_matrix,SR_region_read_count_matrix * SR_region_eff_length_matrix)
        # SR_isoform_expression = SR_expected_counts / isoform_lengths
        LR_isoform_expression = assign_reads(isoform_expression,LR_isoform_region_matrix,LR_region_read_count_matrix)
    prediction_params = None
    if (alpha == 'adaptive'):
        # if gene == 'ENSG00000280987.4':
        #     with open('temp.pkl','wb') as f:
        #         pickle.dump([SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix],f)
        prediction_params = (SR_isoform_region_matrix,SR_region_read_count_matrix,\
            LR_isoform_region_matrix,LR_region_read_count_matrix,\
                isoform_lengths,isoform_num_exons,model)

    # for params in params_grid:

    return LR_isoform_expression,prediction_params
def quantification(short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,SR_gene_isoform_expression_dict,SR_quantification_option,DL_model,alpha,beta,P,assign_unique_mapping_option):
    print('[INFO] Calculating the isoform expression...',flush=True)
    gene_isoform_expression_dict = defaultdict(lambda:defaultdict(dict))
    if (alpha == 'adaptive'):
        model = DL_model
    else:
        model = None
    print(f'[INFO] Using alpha = {alpha}',flush=True)
    list_of_all_genes_chrs = []
    for chr_name in long_read_gene_matrix_dict:
        if chr_name in short_read_gene_matrix_dict:
            for gene_name in long_read_gene_matrix_dict[chr_name]:
                if gene_name in short_read_gene_matrix_dict[chr_name]:
                    list_of_all_genes_chrs.append((gene_name,chr_name))
    list_of_args = [(short_read_gene_matrix_dict[chr_name][gene_name],long_read_gene_matrix_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],alpha,beta,P,model,assign_unique_mapping_option,SR_quantification_option,gene_name) for gene_name,chr_name in list_of_all_genes_chrs]
    print('[INFO] Solving by linear model...',flush=True)
    list_of_prediction_params = []
    for (gene_name,chr_name), args in zip(list_of_all_genes_chrs, list_of_args):
        # try:
        result = estimate_isoform_expression_single_gene(args)
        # if SR_quantification_option == 'Mili':
        #     gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'],\
        #         gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts'],\
        #             gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'],\
        #                 gene_isoform_expression_dict[chr_name][gene_name]['alpha'] = result
        # else:
        gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] = SR_gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression']
            # gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'],\
            #     gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts'] = \
            #         SR_gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'],\
            #         SR_gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts']
        gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'],prediction_params = result
        list_of_prediction_params.append((gene_name,chr_name,prediction_params))
    print('[INFO] Done',flush=True)
    print('[INFO] Predicting alpha by deep learning model...',flush=True)
    list_of_gene_alpha = predict_params_all_genes(list_of_prediction_params)
    for (gene_name,chr_name, alpha) in list_of_gene_alpha:
        gene_isoform_expression_dict[chr_name][gene_name]['alpha'] = alpha
    print('[INFO] Done',flush=True)

    gene_isoform_tpm_expression_dict = normalize_expression(gene_isoform_expression_dict)
    return gene_isoform_tpm_expression_dict,list_of_all_genes_chrs