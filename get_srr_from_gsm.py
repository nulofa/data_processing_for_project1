import numpy as np
s2g = np.load('total_srr2gsm_dict.npy', allow_pickle=True).item()

g2s = {}
for srr in s2g:
    g2s[s2g[srr]] = srr

with open('fake_exp_gsm') as f:
    for line in f:
        if not line.startswith('GSM'):
            continue
        gsm = line.strip()
        print(g2s[gsm])