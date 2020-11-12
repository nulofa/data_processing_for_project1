import numpy as np

#sam_detail2, sam2gsm_final2都是改了key的，改为filename+'---'+colname为key
#data_dels 是删除了部分非样本的数据，dels2删除了无法匹配的部分数据
s2g = np.load('./sam2gsm_final.npy', allow_pickle=True)
s2g = s2g.item()
sam_cond = np.load('./sample_cond_dels2.npy')
sam_d = np.load('./sam_details_dels2.npy', allow_pickle=True)
sam_d = sam_d.item()
data = np.load('./dataset_new_dels2.npy')
i = 0
cnt=0
dels = []
for s in sam_cond:

    if s not in sam_d:
        print(s, i)
        cnt+=1
        dels.append(i)
    i+=1
data = np.delete(data, dels, axis=0)
print(data.shape)
# print(data[0][:10])
print(cnt, i)
sd = np.delete(sam_cond, dels)
print(len(sd))


# np.save('./sample_cond_dels2', sd)
# np.save('./dataset_new_dels2.npy', data)



