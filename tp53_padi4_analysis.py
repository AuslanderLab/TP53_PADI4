# Imports --------------------------------------------------------------------------------------------------------------
import pickle
import pandas as pd
import numpy as np
import lifelines
from lifelines import KaplanMeierFitter
import scipy
from scipy import stats
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.facecolor'] = 'none'

def extract_cancer(expt,clin,canc='SKCM'):
    smp = clin[clin.type == canc].index
    inc = list(set(smp)&set(expt.columns))
    dexp = expt[inc]
    dexp = dexp.loc[:, ~dexp.columns.duplicated()].copy()
    dclin = clin.loc[inc]
    return dexp,dclin

def plot_km(val,surv,death, mkplt = True, THR=0):

    # if val.name is None:
    #     val.name = 'prots'
    THR = np.median(val)
    dat = pd.DataFrame({'surv': list(surv), 'death': list(death), 'val': list(val)})
    lr = lifelines.statistics.logrank_test(dat[dat.val <= THR].surv, dat[dat.val > THR].surv,
                                           event_observed_A=dat[dat.val <= THR].death,
                                           event_observed_B=dat[dat.val > THR].death)

    kmf = KaplanMeierFitter()
    kmf.fit(dat[dat.val <= THR].surv, event_observed=dat[dat.val <= THR].death, label="low, n="+str(len(dat[dat.val <= THR].death)))
    s1 = kmf.median_survival_time_
    if mkplt:


        ax=kmf.plot_survival_function()

    kmf.fit(dat[dat.val > THR].surv, event_observed=dat[dat.val > THR].death, label="high, n="+str(len(dat[dat.val > THR].death)))
    s2 = kmf.median_survival_time_

    if mkplt:
        kmf.plot_survival_function(ax=ax)

        ax.set_ylabel('survival probability')
        ax.set_xlabel('time')
        mytext = "P=%.2e" % (lr.p_value)
        ax.text(0.1, 0.1, mytext)

    return lr.p_value,s1/s2

def get_corr_mat(expt,clin):
    genesc = [ 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']

    R=[];P=[]
    for t in types:
        rr=[];pp=[]
        dexp, dclin = extract_cancer(expt, clin, canc=t)
        for gn in genesc:
            g1=dexp.loc['PADI4']
            g2=dexp.loc[gn]
            r, p = stats.spearmanr(g1,g2)
            rr.append(r)
            pp.append(p)
        R.append(rr)
        P.append(p)
    df=pd.DataFrame(R,columns=genesc,index=types)
    df=df.assign(m=df.mean(axis=1)).sort_values('m').drop('m', axis=1)
    df.to_csv('outs/corr_gene.csv')

def make_scatters(expt,clin):
    plt.clf()
    ##SKCM corr
    dexp, dclin = extract_cancer(expt, clin, canc='SKCM')
    plt1 = plt.scatter(stats.zscore(dexp.loc['PADI4']),
                       stats.zscore(dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean(), c="blue",
                       alpha=0.3)
    plt.xlabel('PADI4')
    plt.ylabel('avg(CEACAM21,IL16,S1PR4,IL21R')
    r, p = scipy.stats.pearsonr(dexp.loc['PADI4'], dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']].mean())
    mytxt = mytext = "R=%.2e;P=%.2e" % (r, p)
    plt.text(-1, -2, mytxt)
    plt.xlim([-0.9, 5.9])
    plt.savefig('outs/skcm_corr.pdf')

    ##COAD corr
    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='COAD')
    plt1 = plt.scatter(stats.zscore(dexp.loc['PADI4']),
                       stats.zscore(dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean(), c="blue",
                       alpha=0.3)
    plt.xlabel('PADI4')
    plt.ylabel('avg(CEACAM21,IL16,S1PR4,IL21R')
    r, p = scipy.stats.pearsonr(dexp.loc['PADI4'], dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']].mean())
    mytxt = mytext = "R=%.2e;P=%.2e" % (r, p)
    plt.text(-1, -2, mytxt)
    plt.xlim([-1.1, 4])
    plt.savefig('outs/coad_corr.pdf')

    ##PAAD corr
    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='PAAD')
    plt1 = plt.scatter(stats.zscore(dexp.loc['PADI4']),
                       stats.zscore(dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean(), c="blue",
                       alpha=0.3)
    plt.xlabel('PADI4')
    plt.ylabel('avg(CEACAM21,IL16,S1PR4,IL21R')
    r, p = scipy.stats.pearsonr(dexp.loc['PADI4'], dexp.loc[['CEACAM21', 'IL16', 'S1PR4', 'IL21R']].mean())
    mytxt = mytext = "R=%.2e;P=%.2e" % (r, p)
    plt.text(-1, -2, mytxt)
    plt.xlim([-1.7, 3.2])
    plt.savefig('outs/paad_corr.pdf')


def plot_km_tcga(expt,clin):
    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='SKCM')
    surv = dclin['DSS.time']
    death = dclin['DSS']
    val = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()
    plot_km(val, surv, death, mkplt=True)
    plt.savefig('outs/skcm_km.pdf')

    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='CESC')
    surv = dclin['DSS.time']
    death = dclin['DSS']
    val = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()
    plot_km(val, surv, death, mkplt=True)
    plt.savefig('outs/skcm_cesc.pdf')

    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='HNSC')
    surv = dclin['OS.time']
    death = dclin['OS']
    val = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()
    plot_km(val, surv, death, mkplt=True)
    plt.savefig('outs/skcm_hnsc.pdf')

def calc_all_corr_cibersort(expt,clin,ciber):
    C = [];
    RC=[]; PC = []
    for t in types:
        dexp, dclin = extract_cancer(expt, clin, canc=t)
        ss = list(set(dexp.columns) & set(ciber.index))
        dexp = dexp[ss]
        sc = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()

        cii = ciber.loc[ss]
        rr = [];
        pp = []
        for c in ciber.columns:
            r, p = stats.spearmanr(sc, cii[c])
            rr.append(r)
            pp.append(p)
        RC.append(rr)
        PC.append(pp)
    df = pd.DataFrame(RC, columns=ciber.columns, index=types)
    df.to_csv('outs/ciberres.csv')

def make_ciber_scatter(expt,clin,ciber):
    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='SKCM')
    ss = list(set(dexp.columns) & set(ciber.index))
    dexp = dexp[ss]
    sc = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()
    cii = ciber.loc[ss]
    plt1 = plt.scatter(cii['T.cells.CD8'], sc, c="red", alpha=0.4)
    r, p = scipy.stats.pearsonr(cii['T.cells.CD8'], sc)
    plt.xlabel('CD8 Tcells')
    plt.ylabel('avg(PADI4,CEACAM21,IL16,S1PR4,IL21R')
    mytxt = "R=%.2e;P=%.2e" % (r, p)
    plt.text(0, -2, mytxt)
    plt.xlim([-0.01, 0.502])
    plt.savefig('outs/cd8_skcm_scatter.pdf')

    plt.clf()
    dexp, dclin = extract_cancer(expt, clin, canc='SKCM')
    ss = list(set(dexp.columns) & set(ciber.index))
    dexp = dexp[ss]
    sc = stats.zscore(dexp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']], axis=1).mean()
    cii = ciber.loc[ss]
    plt1 = plt.scatter(cii['T.cells.regulatory..Tregs.'], sc, c="red", alpha=0.4)
    r, p = scipy.stats.pearsonr(cii['T.cells.regulatory..Tregs.'], sc)
    plt.xlabel('Tregs')
    plt.ylabel('avg(PADI4,CEACAM21,IL16,S1PR4,IL21R')
    mytxt = "R=%.2e;P=%.2e" % (r, p)
    plt.text(0, -2, mytxt)
    plt.xlim([-0.01, 0.2])
    plt.savefig('outs/tregs_skcm_scatter.pdf')


def make_liu_anti_pd1_plot():
    plt.clf()
    exp = pd.read_csv('ici_data/liu_exp.csv', index_col='Hugo_Symbol')
    lc = pd.read_csv('ici_data/liu_clin.csv', index_col='Patient')
    idd = set(lc.index) & set(exp.columns)
    lc = lc.loc[idd]
    exp = exp[idd]
    v = exp.loc[['PADI4', 'CEACAM21', 'IL16', 'S1PR4', 'IL21R']].mean()
    surv = lc['OS']
    dead = lc['dead']
    plot_km(v, surv, dead, mkplt=True)
    plt.savefig('outs/liu_km.pdf')


def make_hugo_anti_pd1_plot():
    plt.clf()
    rexp = pd.read_csv('ici_data/hugo_exp.csv', index_col='Entrez_Gene_Id')
    rc = pd.read_csv('ici_data/hugo_clin.csv', index_col='PATIENT_ID')
    idd = set(rc.index) & set(rexp.columns)
    rc = rc.loc[idd]
    rexp = rexp[idd]
    v = rexp.loc[[23569, 90273, 3603, 8698, 50615]].mean()
    surv = rc['OS_MONTHS']
    dead = rc['OS_STATUS'] == '1:DECEASED'
    plot_km(v, surv, dead, mkplt=True)
    plt.savefig('outs/hugo_km.pdf')


def main():
    #load and process datasets

    ##GE data
    exp = pd.read_csv('tcga_data/tcga_RSEM_Hugo_norm_count.txt',sep='\t',index_col='sample')
    smpt = [i for i in list(exp.columns) if '-11' not in i]
    expt = exp[smpt]
    expt=expt.rename(columns={i:i[:12] for i in list(expt.columns)})
    expt = expt.loc[:, ~expt.columns.duplicated()].copy()

    ##clinical data
    clin = pd.read_csv('tcga_data/tcga_clin.csv',index_col='sample_id')
    inc = list(set(clin.index)&set(expt.columns))
    clin = clin.loc[inc]
    expt = expt[inc]
    clin=clin.fillna(0)

    #cibersort deconv
    ciber=pd.read_csv('tcga_data/tcga_cybersort.csv',index_col='SampleID')
    ciber=ciber.T
    ciber = ciber.loc[:, ~ciber.columns.duplicated()].copy()
    ciber=ciber.T

    types = ['OV', 'CESC', 'UVM', 'BRCA', 'KIRC', 'LGG', 'KICH', 'COAD', 'SKCM', 'ACC', 'ESCA', 'UCS', 'CHOL', 'UCEC', 'THYM', 'READ', 'MESO', 'KIRP', 'SARC', 'STAD', 'PRAD', 'LUSC', 'TGCT', 'LUAD', 'LIHC', 'PCPG', 'GBM', 'HNSC', 'PAAD', 'DLBC', 'THCA', 'BLCA']
    gene_e = [23569,90273,3603,8698,50615]
    gene_s = ['PADI4','CEACAM21','IL16','S1PR4','IL21R']


    ##get data for panel A (which is generated by a seperate R code)
    get_corr_mat(expt,clin)

    ##specific scatters panel B:
    make_scatters(expt,clin)

    ##Make KM TCGA (Panel C)
    plot_km_tcga(expt,clin)


    ##get table with all cibersort associations
    calc_all_corr_cibersort(expt,clin,ciber)

    ## plot cibersort scatters (Panel  D)
    calc_all_corr_cibersort(expt,clin,ciber)

    ## plot anti-PD1 survival plots (Panel E)
    make_liu_anti_pd1_plot()
    make_hugo_anti_pd1_plot()

