import numpy as np

sam2gsm = np.load('./sam2gsm.npy',allow_pickle=True)
sam2gsm0 = np.load('./sam2gsm_manul_1.npy', allow_pickle=True)
sam2gsm0 = sam2gsm0.item()
sam2gsm = sam2gsm.item()
d1=[]
d2=[]
for s in sam2gsm:
    print(s in sam2gsm0)
    d1.append((s, sam2gsm[s]))
print('==================')
for s in sam2gsm0:
    d2.append((s, sam2gsm0[s]))

m = len(d1)
n = len(d2)
print(m,n)
mn = min(m,n)
# for i in range(m):
#     print(d1[i])

for i in range(n):
    print(d1[i],d2[i])