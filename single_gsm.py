import numpy as np

sam_d = np.load('./sam_details_srr+new.npy', allow_pickle=True).item()

cur_gsm = {}
with open('./new_gsm') as f:
    for line in f:
        gsm = line.strip()
        if gsm not in cur_gsm:
            cur_gsm[gsm] = 1

for sc in sam_d:
    ogsm = sam_d[sc][0]
    if ogsm not in cur_gsm:
        cur_gsm[ogsm] = 1

total_srr = []
with open('./srrlist.txt') as f:
    for line in f:
        srr = line.strip()
        total_srr.append(srr)

total_gsm = []
with open('./srr_gsm.txt') as f:
    for line in f:
        gsm = line.strip()
        total_gsm.append(gsm)

dels = []
for i in range(len(total_gsm)):
    if total_gsm[i] in cur_gsm:
        dels.append(i)

se_srr = np.delete(total_srr, dels)
se_gsm = np.delete(total_gsm, dels)

with open('./pe_srr.txt', 'w') as f:
    for srr in se_srr:
        f.write(srr+'\n')

with open('./pe_gsm.txt', 'w') as f:
    for gsm in se_gsm:
        f.write(gsm+'\n')
