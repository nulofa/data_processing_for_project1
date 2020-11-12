import numpy as np
import os

# with open('./new_data0/rsem_res.genes.results') as f:
#     fstr = ''
#     for line in f:
#         cont = line.split('\t')
#         cont[0] = cont[0].replace('gene-', '')
#         fstr += cont[0] + '\t' + cont[-2] + '\n'
#
#     with open('./new_data0/resm.res.txt', 'w') as f:
#         f.write(fstr)

gene_name = np.load('./gene_name2.npy')
gn_dup = np.load('./gn_dup2.npy', allow_pickle=True).item()
gn2ind = {}
gene_len = len(gene_name)
print(gene_len)
for i in range(gene_len):
    gn2ind[gene_name[i]] = i
    if gene_name[i] in gn_dup:
        for gn in gn_dup[gene_name[i]]:
            gn2ind[gn] = i

files = os.listdir('./new_data0/')
new_cond=[]
new_data = []

for fi in files:
    datap = []
    cond = []
    if 'txt' in fi:
        # continue
        with open('./new_data0/'+fi) as f:

            for line in f:
                line = line.rstrip()
                cont = line.split('\t')
                gn = cont[0].strip()
                if gn == '' or gn not in gn2ind:
                    if len(cond) == 0:
                        cond.extend(cont[1:])
                        sample_cnt = len(cont) - 1
                        # print(cont)
                        # print(fi+'\n\n')
                    # print('continue')
                    continue
                # print('no continue')
                if len(datap) == 0:
                    datap = np.zeros((sample_cnt, gene_len + 1), dtype='float32')
                for i in range(sample_cnt):
                    # print(float(cont[i+1].strip()))
                    try:
                        if cont[i + 1].strip() == 'nan':
                            datap[i][gn2ind[gn]] = 0.0
                        else:
                            if float(cont[i + 1].strip()) == np.inf or float(cont[i + 1].strip()) == -np.inf:
                                print('----------------表达量太大')
                                break
                            datap[i][gn2ind[gn]] = float(cont[i + 1].strip())
                    except:
                        print('try error')
        new_data.extend(datap)
        print(fi, np.array(datap).shape)
        for cd in cond:
            if not fi.startswith('GSM'):
                new_cond.append(fi+'---'+cd)
            else:
                new_cond.append(fi.split('_')[0])
    else:
        continue


    new_data.extend(datap)

from scipy.stats.stats import spearmanr
print(spearmanr(new_data[0], new_data[1]))

