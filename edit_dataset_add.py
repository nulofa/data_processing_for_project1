import pickle
import re
import os

import numpy as np
import requests



def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"

sam_d = np.load('./sam_details_dels2.npy', allow_pickle=True).item()
sam_cond = np.load('./sample_cond_dels2.npy')
dataset0 = np.load('./dataset_dels2.npy')
dataset0 = list(dataset0)

data_pq = np.load('./dataset_pq.npy')
gene_name = np.load('./gene_name2.npy')
gn_dup = np.load('./gn_dup2.npy', allow_pickle=True).item()
gn2ind = {}
gene_len = len(gene_name)
print(gene_len)
cnt = 0
for i in range(gene_len):
    gn2ind[gene_name[i]] = i
    if gene_name[i] in gn_dup:
        for gn in gn_dup[gene_name[i]]:
            gn2ind[gn] = i

files = os.listdir('./new_exp/')



cond=[]
new_data = []
fn = []
cnt = 0
for fi in files:
    datap = []
    if 'sample' in fi:
        fn.append(fi.split('_')[1])
    else:
        fn.append(fi.split('.')[0])
    # print(fn)

    with open('./new_exp/'+fi) as f:
        for line in f:
            line = line.rstrip()
            cont = line.split('\t')
            gn = cont[0].strip()

            if gn == '' or gn not in gn2ind:
                if len(cond) == 0:
                    cond.extend(cont[1:])
                    sample_cnt = len(cont) - 1
                    # print(sample_cnt)
                continue
            if len(datap) == 0:
                datap = np.zeros((sample_cnt, gene_len+1), dtype='float32')
            for i in range(sample_cnt):
                # print(float(cont[i+1].strip()))
                try:
                    if cont[i + 1].strip() == 'nan':
                        datap[i][gn2ind[gn]] = 0.0
                    else:
                        if float(cont[i + 1].strip()) == np.inf or float(cont[i + 1].strip()) == -np.inf:
                            print('----------------表达量太大')
                            break
                        datap[i][gn2ind[gn]] = float(cont[i + 1].strip())
                except:
                    print('try error')
    new_data.extend(datap)
    if np.sum(datap) == 0:
        print(fi)


new_n = len(new_data)
ndels = []
i4nd = 0
print(new_n)
for d in np.array(new_data):
    if np.sum(d>1) < 5000:
        ndels.append(i4nd)
    i4nd+=1
new_data = np.delete(new_data, ndels, axis=0)
dataset0.extend(new_data)
np.save('./dataset_10_2',dataset0)
print('len of dataset0 = ',len(dataset0))

fn = np.delete(fn, ndels, axis=0)
print(len(new_data), len(ndels), len(fn))
new_n -= len(ndels)
srr2gsm = np.load('./total_srr2gsm_dict.npy', allow_pickle=True).item()



gsmlist = []
condlist =  []
for srr in fn:
    gsmlist.append(srr2gsm[srr])
    condlist.append(srr)
#
#
#
#
sam_cond = list(sam_cond)
sam_cond.extend(condlist)


np.save('./sample_cond_10_2', sam_cond)
print(sam_cond[-new_n:])
#add to sam_dettail
for sc in sam_cond[-new_n:]:
    print(sc)
    gsm = srr2gsm[sc]
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + gsm
    gsmcont = getHtmlText(url)
    # tissue = re.findall('tissue: ([^:]*)<br>', gsmcont)
    title = re.findall('<tr valign="top"><td nowrap>Title</td>\n<td style="text-align: justify">([^\n]*)</td>\n</tr>',
                       gsmcont)

    # if len(tissue) == 0:
    #     tissue = re.findall('<tr valign="top"><td nowrap>Source name</td>\n<td style="text-align: justify">([^:]*)<br></td>',gsmcont)
    # genotype = re.findall('genotype[a-zA-Z/]*: ([a-zA-Z0-9\s]*)',gsmcont)
    # if len(genotype) == 0:
    #     genotype = re.findall('ecotype: ([^:]*)<br>', gsmcont)
    # if len(genotype) == 0:
    #     genotype = re.findall('genetype: ([^:]*)<br>', gsmcont)
    # treatment = re.findall('<br>treatment: ([^:]*)<br>',gsmcont)
    # if len(treatment) == 0:
    # treatment = re.findall('<tr valign="top"><td nowrap>Treatment protocol</td>\n<td style="text-align: justify">([.]*)<br></td>', gsmcont)
    characteristics = re.findall(
        '<tr valign="top"><td nowrap>Characteristics</td>\n<td style="text-align: justify">(.*)<br></td>', gsmcont)
    # sam_details[gse+ttl] = (gsm, tissue, genotype, treatment)

    # print((gsm, tissue, genotype, treatment))
    sam_d[sc] = (gsm, title, characteristics)
    print((gsm, title, characteristics), sc)

np.save('./sam_details_10_2', sam_d)
print(len(sam_d), len(sam_cond))
def d2rank(data):
    for d in data:
        n0s = np.sum(d<=0)
        d2r = np.argsort(d)
        d[d2r] = range(38313)
        d -= n0s
        d[d<0] = 0
    return data
def standardrized(d):
    mu = np.mean(d)
    sigmoid = np.std(d)

    return (d-mu)/sigmoid

# add to dataset

dataset=np.load('./data_rank_and_norn_10_2.npy')
dataset = list(dataset)
print(len(dataset))
pq = None
with open('pq.pkl', 'rb') as f:
        pq= pickle.load(f)
data_pq = list(data_pq)

new_data = d2rank(new_data)
for d in new_data:
    newd = standardrized(d)
    print(np.sum(newd), np.std(newd))
    dataset.append(newd)
    data_pq.append(pq.encode(np.reshape(newd, (1,-1))))
# np.save('./data_rank_and_norm_10_2', dataset)

print(len(dataset0))
print(len(dataset))
# np.save('./dataset_pq_10_2',data_pq)


# np.save('./current_gse', gselist)



