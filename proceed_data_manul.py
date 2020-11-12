import numpy as np
import os


# new_cond = np.load('./manul_cond.npy')
# new_s2g = np.load('./manul_s2g.npy',allow_pickle=True).item()
#
# for cd in new_cond:
#     print(cd, new_s2g[cd])
#     continue
#     if not cd.startswith('GSM3242'):
#         continue
#     if cd.startswith('GSM'):
#         gsm = cd
#     else:
#         gsm = input('gsm = ')
#     new_s2g[cd] = gsm
#
# # np.save('manul_s2g', new_s2g)
# assert 1==0

# files = os.listdir('./newdata_manul/')
# gene_name = np.load('./gene_name2.npy')
# gn_dup = np.load('./gn_dup2.npy', allow_pickle=True).item()
# gn2ind = {}
# gene_len = len(gene_name)
# for i in range(gene_len):
#     gn2ind[gene_name[i]] = i
#     if gene_name[i] in gn_dup:
#         for gn in gn_dup[gene_name[i]]:
#             gn2ind[gn] = i
#
# new_cond = list(np.load('./manul_cond.npy'))
# new_data = list(np.load('./manul_data.npy'))
# for fi in files:
#     cond = []
#     datap = []
#     # print(fi)
#     with open('./newdata_manul/'+fi) as f:
#         if not fi.startswith('GSM3242'):
#             continue
#         for line in f:
#             line = line.rstrip()
#             cont = line.split('\t')
#             gn = cont[0].strip()
#             if gn == '' or gn not in gn2ind:
#                 if len(cond) == 0:
#                     cond.extend(cont[1:])
#                     sample_cnt = len(cont) - 1
#                     # print(cont)
#                     # print(fi+'\n\n')
#                 # print('continue')
#                 continue
#             # print('no continue')
#             if len(datap) == 0:
#                 datap = np.zeros((sample_cnt, gene_len + 1), dtype='float32')
#             for i in range(sample_cnt):
#                 # print(float(cont[i+1].strip()))
#                 try:
#                     if cont[i + 1].strip() == 'nan':
#                         datap[i][gn2ind[gn]] = 0.0
#                     else:
#                         if float(cont[i + 1].strip()) == np.inf or float(cont[i + 1].strip()) == -np.inf:
#                             print('----------------表达量太大')
#                             break
#                         datap[i][gn2ind[gn]] = float(cont[i + 1].strip())
#                 except:
#                     print('try error')
#     new_data.extend(datap)
#     print(fi, np.array(datap).shape)
#     for cd in cond:
#         if not fi.startswith('GSM'):
#             new_cond.append(fi + '---' + cd)
#         else:
#             new_cond.append(fi.split('_')[0])
#     # print(new_cond)
#     print(len(new_cond), len(new_data))
# np.save('./manul_data',new_data)
# np.save('./manul_cond', new_cond)
# for cd in new_cond:
#     print(cd)

# assert 1==0 # processing original data
gene_length = np.load('./gene_length.npy', allow_pickle=True).item()

def cal_exp(pwd,filename):
    gene_name = np.load(pwd+'/gene_name2.npy')
    gn_dup = np.load(pwd+'/gn_dup2.npy', allow_pickle=True).item()
    gn2ind = {}
    gene_len = len(gene_name)
    print(gene_len)
    for i in range(gene_len):
        gn2ind[gene_name[i]] = i
        if gene_name[i] in gn_dup:
            for gn in gn_dup[gene_name[i]]:
                gn2ind[gn] = i
    cnt_d = []
    sam_n = []
    c4l = 0
    with open(filename) as f:
        for line in f:
            c4l+=1
            # if c4l > 100:
            #     break
            cont = line.split('\t')

            if cont[-1] == '\n':
                cont = cont[:-1]

            if len(sam_n) == 0:
                sam_n.extend(cont[1:])
                print('sam = ', sam_n)
            if len(cont) < 1:
                continue
            gn = cont[0]
            if gn not in gn2ind:
                # print(gn, 'not in gn2ind')
                continue

            cnt = [float(c) for c in cont[1:]]


            if gn not in gene_length:
                gind = gn2ind[gn]
                tmp_gn = []
                for tgn, gi in gn2ind.items():
                    if gi == gind:
                        tmp_gn.append(tgn)
                tmp_gn.sort()
                # print(tmp_gn)
                for tgn in tmp_gn:
                    if tgn.startswith('AT') and len(tgn) > 3 and tgn[3] == 'G':
                        gn = tgn
                        break
            if gn not in gene_length:
                # print(gn, 'not in gene length')
                continue
            cntd = [gn , gene_length[gn]]
            cntd.extend(cnt)
            cnt_d.append(cntd)
            # print(cnt_d)

        counts = []
        for d in cnt_d:

            counts.append(np.array(d[2:]) / np.array(d[1]))

        counts = np.array(counts)
        print(counts.shape, np.sum(counts ,axis=0).shape)
        tpm = counts / np.sum(counts ,axis=0)*1e6

        # count_data = pd.DataFrame(columns=col, data= cnt_d)
        #
        # reads = count_data.loc[:, [sn for sn in sam_n]]
        #
        # genes = count_data.loc[:, 'gene']
        #
        # glength = count_data.loc[:, 'length']
        # rate = np.array(reads.values) / np.array(glength.values)
        # tpm = rate / np.sum(rate, axis=0).reshape(1, -1) * 1e6
        # print(tpm)
        out_path = './newdata_manul/'
        with open(out_path+fi, 'w') as f:
            f.write('gene\t'+ '\t'.join(sam_n)+'\n')
            i4tpm = 0
            for d in cnt_d:
                f.write(d[0]+'\t'+'\t'.join(str(i) for i in tpm[i4tpm]))
                f.write('\n')
                i4tpm += 1


import os
import csv
import xlrd
path = './tmp4GSE111490/'
files = os.listdir(path)
skip = True
for fi in files:
    print(fi)

    # if not fi.startswith('GSM4083452'):
    #     continue
    cal_exp('.', path+fi)

