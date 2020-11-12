import pickle
import re
import os
import xlrd
import numpy as np
import requests
import csv
# list = ['GSE112564_pan_induction_0_6_24hours.txt'] #加上以 GSM37 开头的
#
# for fi in list:
#     new_lines = ''
#     with open('./new_data/' + fi) as f:
#         print(fi)
#         contlist = input('cont list = ').split(' ')
#
#         for line in f:
#             cont = line.split('\t')
#             for ci in contlist:
#                 if ci == contlist[-1]:
#                     new_lines += cont[int(ci)] + '\n'
#                 else:
#                     new_lines += cont[int(ci)] + '\t'
#     with open('./new_data/'+fi, 'w') as f:
#         f.write(new_lines)
# assert 1==0

def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"

sam_cond = np.load('./sample_cond_srr+.npy')
#
sam_d = np.load('./sam_details_srr+.npy', allow_pickle=True).item()
# dataset = np.load('./data_rank_and_norm.npy')
# data_pq = np.load('./dataset_pq.npy')
gene_name = np.load('./gene_name2.npy')
gn_dup = np.load('./gn_dup2.npy', allow_pickle=True).item()
gn2ind = {}
gene_len = len(gene_name)
print(gene_len)
for i in range(gene_len):
    gn2ind[gene_name[i]] = i
    if gene_name[i] in gn_dup:
        for gn in gn_dup[gene_name[i]]:
            gn2ind[gn] = i

files = os.listdir('./new_data/')
new_cond=[]
new_data = []

for fi in files:
    datap = []
    cond = []
    if 'txt' in fi:
        # continue
        with open('./new_data/'+fi) as f:

            for line in f:
                line = line.rstrip()
                cont = line.split('\t')
                gn = cont[0].strip()
                if gn == '' or gn not in gn2ind:
                    if len(cond) == 0:
                        cond.extend(cont[1:])
                        sample_cnt = len(cont) - 1
                        # print(cont)
                        # print(fi+'\n\n')
                    # print('continue')
                    continue
                # print('no continue')
                if len(datap) == 0:
                    datap = np.zeros((sample_cnt, gene_len + 1), dtype='float32')
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
        print(fi, np.array(datap).shape)
        for cd in cond:
            if not fi.startswith('GSM'):
                new_cond.append(fi+'---'+cd)
            else:
                new_cond.append(fi.split('_')[0])
    else:
        # print(fi)
        if not fi.endswith('csv'):
            cont_xlsx = xlrd.open_workbook('./new_data/'+fi)
            sh = cont_xlsx.sheet_by_index(0)
            for ri in range(sh.nrows):
                line = sh.row_values(ri)

                gn = line[0]
                if gn == '' or gn not in gn2ind:
                    if len(cond) == 0:
                        cond = line[1:]
                        sample_cnt = len(cond)
                    continue
                if len(datap) == 0:
                    datap = np.zeros((sample_cnt, gene_len + 1), dtype='float32')
                for i in range(sample_cnt):
                    # print(float(line[i+1].strip()))
                    try:

                        if line[i + 1]== 'nan':
                            datap[i][gn2ind[gn]] = 0.0
                        else:
                            if float(line[i + 1]) == np.inf or float(line[i + 1]) == -np.inf:
                                print('----------------表达量太大')
                                break
                            datap[i][gn2ind[gn]] = float(line[i + 1])
                    except:
                        print('try error xlsx')
        else:
            # continue
            with open('./new_data/'+fi) as f:
                cont_csv = csv.reader(f)
                for line in cont_csv:
                    gn = line[0]
                    if gn == '' or gn not in gn2ind:
                        if len(cond) == 0:
                            cond = line[1:]
                            sample_cnt = len(cond)
                        continue
                    if len(datap) ==0:
                        datap =  np.zeros((sample_cnt, gene_len + 1), dtype='float32')
                    for i in range(sample_cnt):
                        # print(float(cont[i+1].strip()))
                        try:
                            if line[i + 1].strip() == 'nan':
                                datap[i][gn2ind[gn]] = 0.0
                            else:
                                if float(line[i + 1].strip()) == np.inf or float(line[i + 1].strip()) == -np.inf:
                                    print('----------------表达量太大')
                                    break
                                datap[i][gn2ind[gn]] = float(line[i + 1].strip())
                        except:
                            print('try error')

        new_data.extend(datap)
        print(fi, np.array(datap).shape)
        for cd in cond:
            new_cond.append(fi + '---' + cd)

new_n = np.array(new_data).shape[0]
print(new_n, ' <-n->', len(new_cond))
s2g = np.load('./sam2gsm_final.npy',allow_pickle=True).item()
# for d in new_cond:
#     if d in s2g:
#         continue
#     if d.startswith('GSM'):
#         s2g[d] = d
#     else:
#         print(d)
#         s2g[d] = input('gsm = ')
# np.save('./sam2gsm_final', s2g)

sam_cond = list(sam_cond)
sam_cond.extend(new_cond)
sam_cond=np.array(sam_cond)

np.save('./sample_cond_srr+new.npy', sam_cond)

print(sam_cond[-new_n:])
#add to sam_dettail
for sc in sam_cond[-new_n:]:
    print(sc)
    gsm = s2g[sc]
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

np.save('./sam_details_srr+new', sam_d)


# def d2rank(data):
#     for d in data:
#         n0s = np.sum(d<=0)
#         d2r = np.argsort(d)
#         d[d2r] = range(38313)
#         d -= n0s
#         d[d<0] = 0
#     return data
# def standardrized(d):
#     mu = np.mean(d)
#     sigmoid = np.std(d)
#
#     return (d-mu)/sigmoid
#
# ## add to dataset
# dataset = np.load('./data_rank_and_norm.npy')
# data_pq = np.load('./dataset_pq.npy')
# dataset = list(dataset)
# pq = None
# with open('pq.pkl', 'rb') as f:
#         pq= pickle.load(f)
# data_pq = list(data_pq)
#
# new_data = d2rank(new_data)
# for d in new_data:
#     newd = standardrized(d)
#     print(np.sum(newd), np.std(newd))
#     dataset.append(newd)
#     data_pq.append(pq.encode(np.reshape(newd, (1,-1))))
# np.save('./data_rank_and_norm_srr+new', dataset)
# np.save('./dataset_pq_srr+new',data_pq)



