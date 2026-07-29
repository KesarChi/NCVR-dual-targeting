import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
import sys
sys.path.append('/home/shaoqi/BIGCode/')
from brush import statools
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')


file_cov = "/pmaster/shaoqi/data/UKB/noProt_UKB_drug.txt.gz"
root = "/pmaster/shaoqi/proj/GeneHunt/LPA/plot3d/"
file_score_lpa = '/pmaster/shaoqi/proj/GeneHunt/LPA/phewas-clinical/jc-lpa.sscore'
file_score_PCSK9 = '/pmaster/shaoqi/proj/GeneHunt/LPA/phewas-clinical/best.score'
file_valid = '/hwmaster/shaoqi/proj/GeneHunt/LPA/phewas/d118/White-rm2nd.val'
file_white = '/hwmaster/shaoqi/proj/GeneHunt/LPA/phewas/d118/White-rm2nd.eid'
lpa = '/hwmaster/shaoqi/proj/GeneHunt/LPA/data/pheno/lpa.pheno'

cov_list = ['smoking', 'age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
y = ['CAD', 'ldl_a']
other_cov1 = ['simvastatin', 'atorvastatin', 'fluvastatin','pravastatin','rosuvastatin', 'eptastatin']
other_cov2 = ['lipitor_10mg_tablet', 'simvador_10mg_tablet', 'lipostat_10mg_tablet', 'crestor_10mg_tablet', 'zocor_10mg_tablet']


cov = pd.read_csv(file_cov, sep='\t', compression='gzip')
cov = cov.drop_duplicates('eid').reset_index(drop=True)
score_lpa = pd.read_csv(file_score_lpa, sep='\s+')
score_PCSK9 = pd.read_csv(file_score_PCSK9, sep='\s+')
eid_val = pd.read_csv(file_valid, sep='\t', header=None)
eid_white = pd.read_csv(file_white, sep='\t', header=None)
lpa = pd.read_csv(lpa, sep='\t')
# ldlrawo = pd.read_csv(ldlrawo, sep='\t')

# score col rename
score_lpa.columns = ['eid', 'LPA']
score_PCSK9.columns = ['eid', 'PCSK9']
lpa.columns = ['eid', 'lpa']

# merge cov and score
cov = cov[['eid'] + cov_list + y]
df = pd.merge(cov, score_lpa, on='eid', )
df = pd.merge(df, score_PCSK9, on='eid',)
df = pd.merge(df, lpa, on='eid', how='left')

# filter test
test = df[(~df['eid'].isin(eid_val[0])) & (df['eid'].isin(eid_white[0]))].reset_index(drop=True)
del cov, score_lpa, score_PCSK9, df, eid_val, eid_white, lpa

test['LPA-res'] = statools.get_resid(test, 'LPA', cov_list[3:])
test['PCSK9-res'] = statools.get_resid(test, 'PCSK9', cov_list[3:])

# scale
test['LPA-res'] = (test['LPA-res'] - test['LPA-res'].mean()) / test['LPA-res'].std()
test['PCSK9-res'] = (test['PCSK9-res'] - test['PCSK9-res'].mean()) / test['PCSK9-res'].std()

test['ldl_a'] = test['ldl_a'] * 38.67

# version 1: 80*80 grid, 1% step, top 10% as top group, 10%-20% as bot group
# res_beta = []
# res_beta_lpal = []
# res_beta_pcsk9l = []
# res_diff_lpa = []
# res_diff_ldl = []

# for i in range(80):
    
#     top1_lpa = test['LPA-res'].quantile(0.9)
#     bot1_lpa = test['LPA-res'].quantile(0.9 - 0.01*i)
#     drop_lpa = test['LPA-res'].quantile(0.8 - 0.01*i)
#     top1_lpa_rows = test[test['LPA-res'] >= top1_lpa].index
#     bot1_lpa_rows = test[(test['LPA-res'] <= bot1_lpa) & (test['LPA-res'] >= drop_lpa)].index
#     test['top-lpa'] = 0
#     test.loc[top1_lpa_rows, 'top-lpa'] = 1
#     test['bot-lpa'] = 0
#     test.loc[bot1_lpa_rows, 'bot-lpa'] = 1
#     for j in range(80):
#         top1_pcsk9 = test['PCSK9-res'].quantile(0.9)
#         bot1_pcsk9 = test['PCSK9-res'].quantile(0.9 - 0.01*j)
#         drop_pcsk9 = test['PCSK9-res'].quantile(0.8 - 0.01*j)
#         top1_pcsk9_rows = test[test['PCSK9-res'] >= top1_pcsk9].index
#         bot1_pcsk9_rows = test[(test['PCSK9-res'] <= bot1_pcsk9) & (test['PCSK9-res'] >= drop_pcsk9)].index
#         test['top-pcsk9'] = 0
#         test.loc[top1_pcsk9_rows, 'top-pcsk9'] = 1
#         test['bot-pcsk9'] = 0
#         test.loc[bot1_pcsk9_rows, 'bot-pcsk9'] = 1
        
#         test['T-LPA-T-PCSK9'] = test['top-lpa'] * test['top-pcsk9']
#         test['T-LPA-B-PCSK9'] = test['top-lpa'] * test['bot-pcsk9']
#         test['B-LPA-T-PCSK9'] = test['bot-lpa'] * test['top-pcsk9']
#         test['B-LPA-B-PCSK9'] = test['bot-lpa'] * test['bot-pcsk9']
        
#         newtest = test[(test['T-LPA-T-PCSK9']==1) | (test['T-LPA-B-PCSK9']==1) | (test['B-LPA-T-PCSK9']==1) | (test['B-LPA-B-PCSK9']==1)].reset_index(drop=True)
#         res_diff_lpa.append(newtest.loc[newtest['top-lpa']==1, 'lpa'].mean() - newtest.loc[newtest['bot-lpa']==1, 'lpa'].mean())
#         res_diff_ldl.append(newtest.loc[newtest['top-pcsk9']==1, 'ldl_a'].mean() - newtest.loc[newtest['bot-pcsk9']==1, 'ldl_a'].mean())
        
#         model = statools.stat_model(newtest, 'CAD', ['B-LPA-B-PCSK9', 'B-LPA-T-PCSK9', 'T-LPA-B-PCSK9'], cov_list[1:], 'logit')
#         res_beta.append(model.params['B-LPA-B-PCSK9'])
#         res_beta_lpal.append(model.params['B-LPA-T-PCSK9'])
#         res_beta_pcsk9l.append(model.params['T-LPA-B-PCSK9'])
    
# sv = pd.DataFrame({'lpa-diff': res_diff_lpa, 'ldl-diff': res_diff_ldl, 'beta-BB': res_beta, 'beta-lpaB': res_beta_lpal, 'beta-pcsk9B': res_beta_pcsk9l})
# sv.to_csv(root + 'plot3d-3.31.csv', index=False, sep='\t')
    
# version 2
res_lpa_percent = []
res_ldl_percent = []
res_beta = []
res_beta_lpal = []
res_beta_pcsk9l = []
res_diff_lpa = []
res_diff_ldl = []
res_p = []
res_p_lpa = []
res_p_ldl = []
res_se = []
res_se_lpa = []
res_se_ldl = []


for i in range(46):
    top1_lpa = test['LPA-res'].quantile(0.5 + 0.01*i)
    bot1_lpa = test['LPA-res'].quantile(0.5 - 0.01*i)
    top1_lpa_rows = test[test['LPA-res'] >= top1_lpa].index
    bot1_lpa_rows = test[test['LPA-res'] <= bot1_lpa].index
    test['top-lpa'] = 0
    test.loc[top1_lpa_rows, 'top-lpa'] = 1
    test['bot-lpa'] = 0
    test.loc[bot1_lpa_rows, 'bot-lpa'] = 1
    for j in range(46):
        res_lpa_percent.append(0.5 - 0.01*i)
        res_ldl_percent.append(0.5 - 0.01*i)
        top1_pcsk9 = test['PCSK9-res'].quantile(0.5 + 0.01*j)
        bot1_pcsk9 = test['PCSK9-res'].quantile(0.5 - 0.01*j)
        top1_pcsk9_rows = test[test['PCSK9-res'] >= top1_pcsk9].index
        bot1_pcsk9_rows = test[test['PCSK9-res'] <= bot1_pcsk9].index
        test['top-pcsk9'] = 0
        test.loc[top1_pcsk9_rows, 'top-pcsk9'] = 1
        test['bot-pcsk9'] = 0
        test.loc[bot1_pcsk9_rows, 'bot-pcsk9'] = 1
        
        test['T-LPA-T-PCSK9'] = test['top-lpa'] * test['top-pcsk9']
        test['T-LPA-B-PCSK9'] = test['top-lpa'] * test['bot-pcsk9']
        test['B-LPA-T-PCSK9'] = test['bot-lpa'] * test['top-pcsk9']
        test['B-LPA-B-PCSK9'] = test['bot-lpa'] * test['bot-pcsk9']
        newtest = test[(test['T-LPA-T-PCSK9']==1) | (test['T-LPA-B-PCSK9']==1) | (test['B-LPA-T-PCSK9']==1) | (test['B-LPA-B-PCSK9']==1)].reset_index(drop=True)
        res_diff_lpa.append(newtest.loc[newtest['top-lpa']==1,'lpa'].mean() - newtest.loc[newtest['bot-lpa']==1,'lpa'].mean())
        res_diff_ldl.append(newtest.loc[newtest['top-pcsk9']==1, 'ldl_a'].mean() - newtest.loc[newtest['bot-pcsk9']==1, 'ldl_a'].mean())

        # newtest = test[(test['T-LPA-T-PCSK9']==1) | (test['T-LPA-B-PCSK9']==1) | (test['B-LPA-T-PCSK9']==1) | (test['B-LPA-B-PCSK9']==1)].reset_index(drop=True)
        
        model = statools.stat_model(newtest, 'CAD', ['B-LPA-B-PCSK9', 'B-LPA-T-PCSK9', 'T-LPA-B-PCSK9'], cov_list[1:], 'logit')
        res_beta.append(model.params['B-LPA-B-PCSK9'])
        res_beta_lpal.append(model.params['B-LPA-T-PCSK9'])
        res_beta_pcsk9l.append(model.params['T-LPA-B-PCSK9'])
        res_se.append(model.bse['B-LPA-B-PCSK9'])
        res_se_lpa.append(model.bse['B-LPA-T-PCSK9'])
        res_se_ldl.append(model.bse['T-LPA-B-PCSK9'])
        res_p.append(model.pvalues['B-LPA-B-PCSK9'])
        res_p_lpa.append(model.pvalues['B-LPA-T-PCSK9'])
        res_p_ldl.append(model.pvalues['T-LPA-B-PCSK9'])
        
sv = pd.DataFrame({'lpa-percent': res_lpa_percent, 'ldl-percent': res_ldl_percent, 'lpa-diff': res_diff_lpa, 'ldl-diff': res_diff_ldl, 'beta-BB': res_beta, 'beta-lpaB': res_beta_lpal, 'beta-pcsk9B': res_beta_pcsk9l, 'p-BB': res_p, 'p-lpaB': res_p_lpa, 'p-pcsk9B': res_p_ldl, 'se-BB': res_se, 'se-lpaB': res_se_lpa, 'se-pcsk9B': res_se_ldl})
sv.to_csv(root + 'plot3d-5.5.csv', index=False, sep='\t')
