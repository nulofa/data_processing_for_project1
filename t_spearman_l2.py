import numpy as np
from scipy.stats.stats import pearsonr,spearmanr


def l2(v1,v2):
    # return hamming_dist(v1, v2)
    assert len(v1) == len(v2)
    v1 =np.array(v1).astype('float32')
    v2 =np.array(v2)
    return np.sum((v1-v2)**2)**0.5

def d2rank(data,flag=1):
    none_zeros = []
    if flag ==1:
        for d in data:

            n0s = np.sum(d<=0)
            d2r = np.argsort(d)
            d[d2r] = range(1,100+1)
            # d -= n0s
            sum_rank_zeros = np.sum(range(1,n0s+1))

            d[d<=n0s] = sum_rank_zeros / n0s

        return data,none_zeros
    else:
        pass

def standardrized(d):
    mu = np.mean(d)
    sigmoid = np.std(d)
    # print(mu,sigmoid)
    return (d-mu)/sigmoid

np.random.seed(1)
data = np.random.rand(30,100)
data[data<0.5] = 0

res0 = []
res1 = []
i = 0
for d in data:

    i += 1
    if i > 2:

        print(spearmanr(data[i-1],data[i-2])[0])
        res0.append(spearmanr(data[i-1],data[i-2])[0])

data, none_zeros = d2rank(data)
for d in data:
    d[:] = standardrized(d)
print('-----')
i=0
for d in data:
    i += 1
    if i > 2:
        print(1-l2(data[i-1],data[i-2])**2/2/10)
        res1.append(1-l2(data[i-1],data[i-2])**2/2/100)

print(abs(np.array(res0) - res1))
print(np.average(abs(np.array(res0) -res1)))