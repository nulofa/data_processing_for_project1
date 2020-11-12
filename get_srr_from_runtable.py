import numpy as np
gsm2srr = {}
up3 = {}
cnt =0
with open('./runtable') as f:
    for line in f:
        cont = line.split(',')
        srr = cont[0]
        gsm = cont[21]
        Assay_Type = cont[1]
        print(Assay_Type)

        if not Assay_Type.upper().startswith('RNA'):
            continue
        # tissue = cont[-4]
        # genotype = cont[-3]
        # if genotype+tissue not in up3:
        #     up3[genotype+tissue] = 1
        # else:
        #     up3[genotype+tissue] =up3[genotype+tissue]+1
        # if up3[genotype+tissue]>3:
        #     continue
        # print(tissue,'---', genotype)
        if not gsm.startswith('GSM'):
            print(srr, gsm)
            continue
        if gsm in gsm2srr:
            # print(gsm, gsm2srr)
            continue
        gsm2srr[gsm] = srr
print(cnt,len(gsm2srr.keys()))
# np.save('./gsm2srr', gsm2srr)
# gsmlist = []
# with open('./srrlist.txt', 'a') as f:
#     for gsm in gsm2srr:
#         gsmlist.append(gsm)
#         f.write(gsm2srr[gsm]+'\n')
#
# with open('./srr_gsm.txt', 'a') as f:
#     for gsm in gsmlist:
#         f.write(gsm+'\n')

