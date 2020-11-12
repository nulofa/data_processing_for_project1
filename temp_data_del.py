import numpy as np
# dataset = np.load('./data_rank_and_norm_srr+new.npy')
s2d = np.load('./sam_details_srr+new.npy', allow_pickle=True).item()
sam_cond = np.load('./sample_cond_srr+new.npy')
# print(dataset.shape)
print(len(sam_cond))
i = 0
dels = []
for d in sam_cond:
    # if 'bisu' in s2d[d][1][0].lower() or 'chap' in s2d[d][1][0].lower():
    #     print(s2d[d][1], s2d[d][0])
    #     dels.append(i)

    if 'isof' in d:
        print(d)

    i+=1

# data_pq = np.load('./dataset_pq_srr+new.npy')
# data_pq = np.delete(data_pq, dels, axis=0)
# np.save('./dataset_pq_srr+new', data_pq)
# sam_cond = np.delete(sam_cond, dels)
# np.save('./sample_cond_srr+new', sam_cond)
# dataset = np.delete(dataset, dels, axis=0)
# np.save('./data_rank_and_norm_srr+new', dataset)
# print(dataset.shape)
# print(len(sam_cond))

# n = dataset.shape[0]
#
# for i in range(10):
#     j = np.random.randint(n)
#     print(dataset[j][:5])
#     print(sam_cond[j])
#     print(s2d[sam_cond[j]])