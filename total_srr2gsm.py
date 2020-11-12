import numpy as np
srr =[]
with open('./total_srr', encoding='utf-8') as f:
    for line in f:
        if not line.startswith('SRR'):
            continue
        srr.append(line.strip())
gsm = []
with open('./totol_gsm', encoding='utf-8') as f:
    for line in f:
        if not line.startswith('GSM'):
            continue
        gsm.append(line.strip())
print(len(srr), len(gsm))

srr2gsm = {}
for i in range(len(srr)):
    srr2gsm[srr[i]] = gsm[i]

np.save('./total_srr2gsm_dict', srr2gsm)

