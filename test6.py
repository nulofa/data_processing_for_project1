import numpy as np
import matplotlib.pyplot as plt
sam_d = np.load('./sam_details_treatment.npy', allow_pickle=True)
sam_d = sam_d.item()
res_list = []
tissues = []
genotypes = []
treatments = []
for s in sam_d:
    gsm, tissue, genotype, treatment = sam_d[s]
    if len(tissue)>0 and len(genotype)>0 and len(treatment)>0:
        res_list.append((s,tissue, genotype, treatment))
        tissues.append(tissue[0])
        genotypes.append(genotype[0])
        treatments.append(treatment[0])
        print(res_list[-1])

tiss2int = {}
i = 0

for ind in np.argsort(tissues):
    tiss2int[tissues[ind]] = i
    i+=1

gt2int = {}
i = 0

for ind in np.argsort(genotypes):
    gt2int[genotypes[ind]] = i
    i+=1

tm2int = {}
i = 0
for ind in np.argsort(treatments):
    tm2int[treatments[ind]] = i
    i+=1

from mpl_toolkits.mplot3d import Axes3D

fig, ax = plt.subplots()
ax = Axes3D(fig)
topk = ['GSE101567ColandINA-RPKM','GSE101567ColandINA-RPKM','GSE101567ColandINA-RPKM','GSE101567Ht1-5andINA-RPKM','GSE101567Nr1andINA-RPKM','GSE101567Col-RPKM']


for r in res_list:
    s, tissue, genotype, treatment = r
    x = tiss2int[tissue[0]]
    y = gt2int[genotype[0]]
    z = tm2int[treatment[0]]
    if s in topk:
        print(s)
        ax.scatter3D(x, y, z, c='r', s=10)
    else:
        ax.scatter3D(x,y,z, c='b', s=10)

plt.savefig('./3d.jpg', dpi=300,bbox_inches='tight')




