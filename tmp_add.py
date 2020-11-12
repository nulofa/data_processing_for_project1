import numpy as np
import os
files = os.listdir('./tmpdata/')

gene_name = np.load('./gene_name2.npy')
gn_dup = np.load('./gn_dup2.npy', allow_pickle=True).item()
gn2ind = {}
gene_len = len(gene_name)
print(gene_len)
cnt = 0
for i in range(gene_len):
    gn2ind[gene_name[i]] = i
    if gene_name[i] in gn_dup:
        for gn in gn_dup[gene_name[i]]:
            gn2ind[gn] = i
cond = []
new_data = []
for fi in files:
    datap = []


    with open('./tmpdata/'+fi) as f:
        for line in f:
            line = line.rstrip()
            cont = line.split('\t')
            gn = cont[0].strip()

            if gn == '' or gn not in gn2ind:
                if len(cond) == 0:
                    cond.extend(cont[1:])
                    sample_cnt = len(cont) - 1
                    # print(sample_cnt)
                continue
            if len(datap) == 0:
                datap = np.zeros((sample_cnt, gene_len+1), dtype='float32')
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
    if np.sum(datap) == 0:
        print(fi)


ndels = []
i4nd = 0
for d in np.array(new_data):
    if np.sum(d>1) < 5000:
        ndels.append(i4nd)
    i4nd+=1
new_data = np.delete(new_data, ndels, axis=0)
dataset = np.load('./dataset_new_dels2.npy')
dataset = list(dataset)
dataset.extend(new_data)
print(len(dataset))
np.save('dataset_new_dels2_2', dataset)