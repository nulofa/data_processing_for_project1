# import nanopq
import numpy as np
import pickle
import re
sam_c = np.load('./sample_cond_srr+new.npy')
sam_d = np.load('./sam_details_srr+new.npy', allow_pickle=True)
sam_d = sam_d.item()
with open('pq.pkl', 'rb') as f:
    pq= pickle.load(f)  # pq_dumped is identical to pq

dataset0 = np.load('./data_rank_and_norm_srr+new.npy')
print(dataset0.shape)
dataset = np.load("./dataset_pq_srr+new.npy")

n = dataset.shape[0]

def edit_distance(str1, str2):
    matrix = [[i + j for j in range(len(str2) + 1)] for i in range(len(str1) + 1)]

    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            if str1[i - 1] == str2[j - 1]:
                d = 0
            else:
                d = 1
            matrix[i][j] = min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1, matrix[i - 1][j - 1] + d)

    return matrix[len(str1)][len(str2)]

#similar
for i in range(n):

    dists = pq.dtable(query=dataset0[i]).adist(codes=dataset.reshape(n,-1))
    min_n = np.argsort(dists)

    gsm = sam_d[sam_c[i]][0]
    # print(gsm)
    for j in min_n[:10]:
        sim = 1 - dists[j] / 38313 / 2
        gsm2 = sam_d[sam_c[j]][0]
        if gsm[:6] == gsm2[:6] or sim <0.98:
            # print(gsm, gsm2)
            continue
        # print(i, sam_d[sam_c[i]])
        print(sam_d[sam_c[i]], sim ,sam_d[sam_c[j]])

    print('\n')
#dis_sim
# for i in range(n):
#
#     dists = pq.dtable(query=dataset0[i]).adist(codes=dataset.reshape(n,-1))
#     min_n = np.argsort(dists)
#     if 'tissue' not in sam_d[sam_c[i]][2][0]:
#         continue
#     else:
#         charac = sam_d[sam_c[i]][2][0].split('<br>')
#         for c in charac:
#             if 'tissue' in c:
#                 tstr = c
#                 break
#     tissue = re.findall('tissue:(.*)',tstr)
#
#     if len(tissue) < 1:
#         print(tstr)
#         break
#     # print(gsm)
#     for j in min_n[-20:-10]:
#         sim = 1 - dists[j] / 38313 / 2
#         if 'tissue' not in sam_d[sam_c[j]][2][0]:
#             continue
#         else:
#             charac2 = sam_d[sam_c[j]][2][0].split('<br>')
#             for c in charac2:
#                 if 'tissue' in c:
#                     tstr2 = c
#                     break
#         tissue2 = re.findall('tissue:(.*)',tstr2)
#         # print(tissue, tissue2)
#         # continue
#
#         if not (tissue[0][:-1] in tissue2[0] or tissue2[0][:-1] in tissue[0]):
#             # print(tissue2, tissue, tissue[0][:-1] in tissue2[0] or tissue2[0][:-1] in tissue[0])
#             continue
#
#         if sim >0.2:
#             # print(gsm, gsm2)
#             continue
#         # print(i, sam_d[sam_c[i]])
#         print(sam_d[sam_c[i]], sim ,sam_d[sam_c[j]])
#
#
#     print('\n')