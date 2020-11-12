import numpy as np
sam2gsm = np.load('./sam2gsm.npy', allow_pickle=True)
sam2gsm = sam2gsm.item()
sam_c = np.load('./sample_cond_new.npy')

i = 0
for c in sam_c:
    gse = c.split('---')[0].split('_')[0]
    col = c.split('---')[1]
    if gse+'---'+col in sam2gsm:
        i +=1
    else:
        print(gse+'---'+col)
        gsm_in = input("gsm = ")
        if gsm_in == 'n':
            continue
        else:
            sam2gsm[gse+'---'+col] = gsm_in
np.save('./sam2gsm_manul.npy', sam2gsm)
