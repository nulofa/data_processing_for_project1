import numpy as np

# sam2gsm = np.load('./sam2gsm_manul_1.npy',allow_pickle=True)
# sam2gsm = sam2gsm.item()
#
# g2s = {}
# dels = []
# for s in sam2gsm:
#     g = sam2gsm[s]
#     if g in g2s:
#         print(s,'----', g, '---',g2s[g])
#         new_g = input('new gsm = ')
#         if new_g == 'n':
#             dels.append(s)
#         else:
#             sam2gsm[s] = new_g
#     else:
#         g2s[g] = s
# print(dels)
# for d in dels:
#     sam2gsm.pop(d)
# np.save('./sam2gsm_final.npy', sam2gsm)

# sam2g = np.load('./sam2gsm_final.npy',allow_pickle=True)
# sam2g = sam2g.item()
# sams = ['GSE115951---Sample_18_PAP_ATP','GSE97857---FPKM.3D1','GSE97857---FPKM.3D2','GSE97857---FPKM.4D1','GSE97857---FPKM.4D2'
#         ,'GSE48661---Col_E_2','GSE48661---Col_AE_2','GSE79576---Col-H']
# gsms = ['GSM3195946','GSM2579256','GSM2579257','GSM2579258','GSM2579259','GSM1183036','GSM1183037','GSM2098530']
# i=0
# for s in sams:
#     sam2g[s] = gsms[i]
#     i+=1
# np.save('./sam2gsm_final1.npy', sam2g)

sam_d = np.load('./sam_details_treatment.npy', allow_pickle=True)
sam_d = sam_d.item()
ip=0
im=0
jiawen = []
mine = []
for s in sam_d:
    gsm, tissue,genotype, treatment = sam_d[s]
    a=len(tissue)==0
    b=len(genotype)==0
    c=len(treatment)==0
    a = int(a)
    b = int(b)
    c = int(c)
    if b==1:
        # print(gsm,'+++++',sam_d[s])
        jiawen.append(gsm)
        ip+=1
    elif a+b+c ==0:
        # print(sam_d[s],'!!!!!!')
        continue
    else:
        print(gsm, '----',sam_d[s])
        im+=1
        mine.append(gsm)
print(ip, im)
#
for g in jiawen:
    print(g)
print('--------------')
for i in mine:
    print(i)

