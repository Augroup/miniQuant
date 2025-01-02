from machine_learning.extract_features import extract_features
from machine_learning.concat_features import concat_features
import config
import pickle
import pandas as pd
import time
import datetime
import numpy as np
import glob
def predict_alpha(output_path,num_SRs,num_LRs):
    all_features_label = ['Num_genes','Num_isoforms','Num_exons','Num_regions',\
  'Min_gene_length','Max_gene_length','Mean_gene_length','Variance_gene_length','Skewness_gene_length','Kurtosis_gene_length',\
  'Min_isoform_lengths','Max_isoform_lengths','Mean_isoform_lengths','Variance_isoform_lengths','Skewness_isoform_lengths','Kurtosis_isoform_lengths',\
  'Min_exon_lengths','Max_exon_lengths','Mean_exon_lengths','Variance_exon_lengths','Skewness_exon_lengths','Kurtosis_exon_lengths',\
 'Min_region_lengths','Max_region_lengths','Mean_region_lengths','Variance_region_lengths','Skewness_region_lengths','Kurtosis_region_lengths',\
  'Num_unique_mapping_SRs','Num_multi_mapping_SRs','Num_mapping_SRs','Num_unique_mapping_LRs','Num_multi_mapping_LRs','Num_mapping_LRs',\
  'Min_unique_mapping_SRs','Max_unique_mapping_SRs','Mean_unique_mapping_SRs','Variance_unique_mapping_SRs','Skewness_unique_mapping_SRs','Kurtosis_unique_mapping_SRs',\
  'Min_multi_mapping_SRs','Max_multi_mapping_SRs','Mean_multi_mapping_SRs','Variance_multi_mapping_SRs','Skewness_multi_mapping_SRs','Kurtosis_multi_mapping_SRs',\
  'Min_unique_mapping_LRs','Max_unique_mapping_LRs','Mean_unique_mapping_LRs','Variance_unique_mapping_LRs','Skewness_unique_mapping_LRs','Kurtosis_unique_mapping_LRs',\
  'Min_multi_mapping_LRs','Max_multi_mapping_LRs','Mean_multi_mapping_LRs','Variance_multi_mapping_LRs','Skewness_multi_mapping_LRs','Kurtosis_multi_mapping_LRs',\
  'Min_5_end_truncation','Max_5_end_truncation','Mean_5_end_truncation','Variance_5_end_truncation','Skewness_5_end_truncation','Kurtosis_5_end_truncation',\
  'Min_3_end_truncation','Max_3_end_truncation','Mean_3_end_truncation','Variance_3_end_truncation','Skewness_3_end_truncation','Kurtosis_3_end_truncation',\
  'Min_SR_kvalue','Max_SR_kvalue','Mean_SR_kvalue','Variance_SR_kvalue','Skewness_SR_kvalue','Kurtosis_SR_kvalue']
    selected_features_index = []
    for l in all_features_label:
        if (('Variance' not  in l) and ('Skewness' not in l) and ('Kurtosis' not in l)):
            selected_features_index.append(all_features_label.index(l))
    print('[INFO] Extracting features to predict the alpha...',flush=True)
    start_time = time.time()
    extract_features(output_path)
    all_features_arr,all_community_ids = concat_features(output_path)
    SR_diff = []
    LR_diff = []
    model_list = []
    for model_path in glob.glob(config.pretrained_model_path+'/*.pkl'):
        coverage = model_path.split('/')[-1]
        LR_cov = int(coverage.split("LR")[-1].split('M_')[0])*10**6
        SR_cov = int(coverage.split("SR")[-1].split('M.')[0])*10**6
        SR_diff.append(num_SRs - SR_cov)
        LR_diff.append(num_LRs - LR_cov)
        model_list.append(model_path)
    cov_diff = np.abs(np.array([SR_diff,LR_diff]).T)
    close_LR_diff = cov_diff[cov_diff[:,1] == np.min(cov_diff[:,1])]
    selected_model = model_list[np.where((cov_diff[:,0] == np.min(close_LR_diff[:,0])) & (cov_diff[:,1]==close_LR_diff[:,1][0]))[0][0]]
    print(f'[INFO] Num_SR:{num_SRs}')
    print(f'[INFO] Num_LR:{num_LRs}')
    print('[INFO] Using pretrained model in {} to predict the alpha...'.format(selected_model),flush=True)
    with open(selected_model,'rb') as f:
        [imputer,clf] = pickle.load(f)
    X = imputer.transform(all_features_arr)
    X = all_features_arr[:,selected_features_index]
    y_pred = clf.predict(X)
    alpha = np.array([i/20 for i in range(0,21)])
    true_strategy_df = pd.DataFrame({'id':all_community_ids,'error':np.min(y_pred,axis=1),'alpha':alpha[np.argmin(y_pred,axis=1)]})
    true_strategy_df.to_csv(config.alpha_df_path,sep='\t',index=False)
    end_time = time.time()
    print('[INFO] Done in {} seconds at {}'.format(end_time-start_time,str(datetime.datetime.now())),flush=True)
    