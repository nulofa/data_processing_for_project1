import struct

import numpy as np

s2g = np.load('./sam2gsm_final.npy', allow_pickle=True)
s2g = s2g.item()
sam_cond = np.load('./sample_cond_dels.npy')
sinf = np.load('./sample_inf.npy')
dataset = np.load('dataset_new_dels.npy')

# di =0
# dels = []
# for d in dataset:
#     if np.sum(d) == 0:
#         dels.append(di)
#         print(sam_cond[di])
#     di+=1
# sam_cond = np.delete(sam_cond, dels)
# np.save('./sample_cond_dels', sam_cond)
# dataset = np.delete(dataset, dels, axis=0)
# np.save('./dataset_new_dels', dataset)
# assert 1==0

dels = []
i=0
for sc in sam_cond:
    # if sc.split('_')[0] == 'GSE74972':
    #     print(sc, i)
    #     dels.append(i)
    # if sc.split('_')[0]+'---'+ sc.split('---')[1]== 'GSE41766---baseMean':
    #     print(sc, i)
    #     dels.append(i)
    # if sc.split('_')[0]+'---'+ sc.split('---')[1] == 'GSE56811---logCPM':
    #     print(sc, i)
    #     dels.append(i)
    if sc.startswith('GSE102713_BPC'):
        print(sc)
        dels.append(i)
    if sc.startswith('GSE56811_full_table_Col0'):
        dels.append(i)
        print(sc)
    if sc == 'GSE80565_RNAseq_ABAtimeSeries_rpkm_log2FC_FDR_filterCPM2n2---AvsE01h_logFC':
        for k in range(i, i+14):
            dels.append(k)
            print(sam_cond[k])
    if sc == 'GSE80567_RNAseq_Dex04h_rpkm_log2FC_FDR_filterCPM2n2---Dex4hDTAF1_logFC':
        for k in range(i, i+4):
            dels.append(k)
            print(sam_cond[k])
    if sc == 'GSE80567_RNAseq_Dex10d_rpkm_log2FC_FDR_filterCPM2n2---Dex10dDTAF1_logFC':
        for k in range(i, i+2):
            dels.append(k)
            print(sam_cond[k])
    if sc.startswith('GSE85847_matrix_rnaseq_GEO_submission_jsantos20160817'):
        dels.append(i)
        print(sc)
    if sc.startswith('GSE67322_RNA_seq'):
        dels.append(i)
        print(sc)
    if sc.startswith('GSE79523') and sc.endswith('start') or sc.endswith('end'):
        dels.append(i)
        print(sc)
    i+=1
print(len(dels))

# sam_cond = np.delete(sam_cond, dels)
# np.save('./sample_cond_dels', sam_cond)
# dataset = np.delete(dataset, dels, axis=0)
# np.save('./dataset_new_dels', dataset)
# print(dataset[j])
# dataset=np.delete(dataset, j, axis=0)

