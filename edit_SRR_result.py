import os
import numpy as np
files = os.listdir('./exp')

for fi in files:
    fstr = ''
    with open('./exp/' + fi) as f:
        for line in f:
            cont = line.split('\t')
            cont[0]= cont[0].replace('gene-', '')
            fstr += cont[0] + '\t' + cont[-2] + '\n'
    with open('./new_exp/'+ fi, 'w') as f:
        f.write(fstr)

# srrlist = []
# with open('./srrlist.txt') as f:
#     for line in f:
#         if not line.startswith('SRR'):
#             continue
#         srrlist.append(line.strip())
# gsmlist = []
# with open('./srr_gsm.txt') as f:
#     for line in f:
#         if not line.startswith('GSM'):
#             continue
#         gsmlist.append(line.strip())
# srr2gsm = np.load("./srr2gsm.npy", allow_pickle=True).item()
# for i in range(len(gsmlist)):
#     srr2gsm[srrlist[i]] = gsmlist[i]
# np.save('./srr2gsm', srr2gsm)
# for fn in files:
#     srr = fn.split('.')[0]
#     print(srr)
