import numpy as np

sam2gsm = np.load('./sam2gsm_manul.npy',allow_pickle=True)
sam2gsm = sam2gsm.item()
sam = ['GSE118102---WOX5_1','GSE118102---WOX5_2','GSE118102---WOX5_3','GSE122355---CR-N50-RPKM',
       'GSE122355---7R-MS-RPKM','GSE122355---7R-N50-RPKM','GSE124340---AWW-1','GSE124340---AWW-2','GSE124340---ADT2-1',
       'GSE124340---ADT2-2','GSE124340---ADT3-1','GSE124340---ADT3-2','GSE124340---ADT4-1','GSE124340---ADT4-2',
       'GSE124340---ADT5-1','GSE124340---ADT5-2','GSE126624---WT_V_1','GSE126624---WT_V_2','GSE126624---WT_V_3',
       'GSE136635---CS-MS-RPKM','GSE48661---Col_AE-1','GSE61545---clf.siliques.r1_FPKM','GSE61545---clf.siliques.r2_FPKM',
       'GSE61545---clf.siliques.r3_FPKM','GSE64381---QC_celseq_GCTTGTC','GSE64381---QC_celseq_GTCGTAT','GSE64381---QC_celseq_GTCGTCG',
       'GSE64381---QC_celseq_GTCGTGA','GSE64381---QC_celseq_GTCGTTC','GSE64381---QC_celseq_TAGGTCG','GSE64381---QC_celseq_TAGGTGA',
       'GSE64381---QC_celseq_TAGGTTC','GSE64381---QC_celseq_TCACACG','GSE64381---QC_celseq_TCACAGA','GSE64381---QC_celseq_TCACATC',
       'GSE64381---QC_celseq_TGATGAT','GSE64381---QC_celseq_TGATGCG','GSE64381---QC_celseq_TGATGGA','GSE64381---QC_celseq_TGATGTC',
       'GSE77211---wt_RPKM','GSE77211---ago1-27_RPKM','GSE80188---CK_FPKM','GSE81202---"AveExp_C.W.0"','GSE81202---"AveExp_C.FR.15"',
       'GSE94995---PFD700Rep1Lane1','GSE94995---PFD100Rep2Lane1','GSE94995---PFD100Rep2Lane2','GSE94995---PFD250Rep1Lane2',
       'GSE98097---fpkm PUB25 Proximal','GSE121673---hos15-r2_hos15-r3','GSE123620---At.col0.u.12.3','GSE123620---At.col0.n.00.1',
       'GSE142537---T3R2p']
gsm = ['GSM3318609','GSM3318610','GSM3318611', 'GSM3464520','GSM3464521','GSM3464522','GSM3529894','GSM3529895','GSM3529896',
       'GSM3529897','GSM3529898','GSM3529899','GSM3529900','GSM3529901','GSM3529902','GSM3529903','GSM3609965','GSM3609966','GSM3609967',
       'GSM4053737','GSM1183037','GSM1507965','GSM1507966','GSM1507967','GSM1570001','GSM1570002','GSM1570003','GSM1570004'
       ,'GSM1570005','GSM1570006','GSM1570007','GSM1570008','GSM1570009','GSM1570010','GSM1570011','GSM1570012','GSM1570013'
       ,'GSM1570014','GSM1570015','GSM2046083','GSM2046084','GSM2114480','GSM2144569','GSM2144577','GSM2494107','GSM2494097',
       'GSM2494098','GSM2494100','GSM2587095','GSM3443227','GSM3508154','GSM3508146','GSM4231551']

dels = ['GSE93420----DSK2 RNAi_control_avg','GSE41766---baseMean','GSE56811---logCPM']
print(len(sam), len(gsm))
for d in dels:
    print(sam2gsm.pop(d))
for i in range(len(gsm)):
    sam2gsm[sam[i]]=gsm[i]
# np.save('./sam2gsm1.npy', sam2gsm)