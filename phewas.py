'''
Author: Kesar-Chi && momiji0613@163.com
Date: 2025-01-14 15:56:33
LastEditors: Kesar-Chi && momiji0613@163.com
LastEditTime: 2025-06-03 15:27:37
FilePath: /BIGCode/project/GeneHunt/LPA/report1/phewas.py
Description: Institute: Beijing Institute Of Genomics, CAS.

Copyright (c) 2025 by ${momiji0613@163.com}, All Rights Reserved. 
'''


import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
from glob import glob
import argparse
import os
import sys
from brush import statools
import warnings
warnings.filterwarnings('ignore')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='disease118')
    parser.add_argument('-lpa', '--lpa', action='store_true', help='lpa')
    parser.add_argument('-pcsk9', '--pcsk9', action='store_true', help='pcsk9')
    args = parser.parse_args()
    
    file_cov = "~/data/noProt_UKB_drug.txt"
    root = "~/proj/phewas/split/"
    if args.lpa:
        save_dir = f"{root}/lpa"
    elif args.pcsk9:
        save_dir = f"{root}/ldl"
    else:
        save_dir = f"{root}/joint"
    file_score_lpa = '~/proj/FIN/LPA/cover0702/model_res/best.score'
    file_score_PCSK9 = '~/proj/decode/cover2/model_res/best.score'
    # file_score_lpa = root + '/jc-lpa.sscore'
    # file_score_PCSK9 = root + '/best.score'
    file_valid = '~/proj/phewas/d118/White-rm2nd.val'
    file_white = '~/proj/phewas/d118/White-rm2nd.eid'
    file_pheno = '~/data/phewas/'

    cov_list = ['smoking', 'age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
    y = ['survive_time', 'incident', 'prevalence', 'CAD', 'tc_a', 'ldl_a', 'apob', 'hdl_a', 'tg_a']
    other_cov1 = ['simvastatin', 'atorvastatin', 'fluvastatin','pravastatin','rosuvastatin', 'eptastatin']
    other_cov2 = ['lipitor_10mg_tablet', 'simvador_10mg_tablet', 'lipostat_10mg_tablet', 'crestor_10mg_tablet', 'zocor_10mg_tablet']


    cov = pd.read_csv(file_cov, sep='\t')
    cov = cov.drop_duplicates('eid').reset_index(drop=True)
    score_lpa = pd.read_csv(file_score_lpa, sep='\s+')
    score_PCSK9 = pd.read_csv(file_score_PCSK9, sep='\s+')
    eid_val = pd.read_csv(file_valid, sep='\t', header=None)
    eid_white = pd.read_csv(file_white, sep='\t', header=None)
    phenos = glob(file_pheno + '*.feather')

    assert score_lpa.shape[1] == 2
    assert score_PCSK9.shape[1] == 2
    score_lpa.columns = ['eid', 'LPA']
    score_PCSK9.columns = ['eid', 'PCSK9']

    # merge cov and score
    cov = cov[['eid'] + cov_list + y + other_cov1 + ['BMI']]
    df = pd.merge(cov, score_lpa, on='eid')
    df = pd.merge(df, score_PCSK9, on='eid')
    print(df.shape)

    test = df[(~df['eid'].isin(eid_val[0])) & (df['eid'].isin(eid_white[0]))].reset_index(drop=True)
    print(test.shape)

    del cov, score_lpa, score_PCSK9, eid_val, eid_white, df
    
    test['LPA'] = statools.get_resid(test, 'LPA', cov_list[3:])
    test['PCSK9'] = statools.get_resid(test, 'PCSK9', cov_list[3:])

    # scale
    test['LPA'] = (test['LPA'] - test['LPA'].mean()) / test['LPA'].std()
    test['PCSK9'] = (test['PCSK9'] - test['PCSK9'].mean()) / test['PCSK9'].std()

    test['LPA-q'] = statools.get_quartile(test, 'LPA', bin=2)
    test['PCSK9-q'] = statools.get_quartile(test, 'PCSK9', bin=2)

    test = statools.get_dummy(test, 'LPA-q')
    test = statools.get_dummy(test, 'PCSK9-q')

    test['LPA-low-PCSK9-low'] = test['LPA-q_1'] * test['PCSK9-q_1']
    test['LPA-low-PCSK9-high'] = test['LPA-q_1'] * test['PCSK9-q_2']
    test['LPA-high-PCSK9-low'] = test['LPA-q_2'] * test['PCSK9-q_1']
    test['LPA-high-PCSK9-high'] = test['LPA-q_2'] * test['PCSK9-q_2']
    
    x = []
    if args.lpa:
        x = ['LPA-q_1']
    elif args.pcsk9:
        x = ['PCSK9-q_1']
    else:
        x = ['LPA-low-PCSK9-low', 'LPA-low-PCSK9-high', 'LPA-high-PCSK9-low']
    
    res = pd.DataFrame({'Trait': [], 'beta': [], 'se': [], 'p': []})
    for i in phenos:
        try:
            ph = pd.read_feather(i)
            ph = ph[['eid', 'event']]
            ph['event'] = ph['event'].astype(float)
            ph['eid'] = ph['eid'].astype(int)
            tmp = pd.merge(test, ph, on='eid')
            if len(tmp['event'].unique()) == 1:
                print(f"{i} has only one unique value")
                
            res_tmp = statools.stat_model(tmp, 'event', x, cov_list[1:], 'logit', verbose=False, save=f"{save_dir}/{i.split('/')[-1].split('.')[0]}.model")
            res = pd.concat([res, pd.DataFrame({'Trait': [i.split('/')[-1].split('.')[0]], 'beta': [res_tmp.params[x[0]]], 'se': [res_tmp.bse[x[0]]], 'p': [res_tmp.pvalues[x[0]]]})])
            if res_tmp.pvalues[x[0]] < 0.05:
                print(f"{i}, {res_tmp.params[x[0]]}, {res_tmp.bse[x[0]]}, {res_tmp.pvalues[x[0]]}")
                
            del tmp, ph
        except Exception as e:
            print(f"{i} has error: {e}")
            continue
        
    res.to_csv(f"{save_dir}/res-phewas.csv", index=False, sep='\t')
        
