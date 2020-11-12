import numpy as np
import os
import csv
import requests
# process xls file
import xlrd
import re

new_sam_d = np.load('./manul_sam_d_total.npy', allow_pickle=True).item()
new_data = np.load('./manul_data_total.npy')
new_cond = np.load('./manul_cond_total.npy')
data = np.load('dataset_10_2.npy')
cond = np.load('sample_cond_10_2.npy')
sam_d = np.load('./sam_details_10_2.npy', allow_pickle=True).item()
data =list(data)
data.extend(list(new_data))
cond =list(cond)
cond.extend(list(new_cond))
print(len(data), len(cond))

for sc in new_sam_d:
    if sc not in sam_d:
        sam_d[sc] = new_sam_d[sc]
    else:
        print(sc)

np.save('dataset_10_21', data)
np.save('sam_cond_10_21', cond)
np.save('sam_details_10_21', sam_d)
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
rankdata = d2rank(data)
std_data = []
for d in rankdata:
    newd = standardrized(d)
    std_data.append(newd)
print(len(std_data))
np.save('./data_rank_and_std_10_21', std_data)
assert 1==0
def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"

sam_d = {}

s2g = np.load('./manul_s2g_total.npy', allow_pickle=True).item()
for sc in s2g:
    if not s2g[sc].startswith('GSM'):
        print(sc, s2g[sc])
    else:
        gsm = s2g[sc]
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + gsm
        gsmcont = getHtmlText(url)
        # tissue = re.findall('tissue: ([^:]*)<br>', gsmcont)
        title = re.findall(
            '<tr valign="top"><td nowrap>Title</td>\n<td style="text-align: justify">([^\n]*)</td>\n</tr>',
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
# np.save('manul_sam_d_total', sam_d)
assert 1==0
#process altered data
files = os.listdir('./new_norms')
files.sort()

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

new_cond = np.load('./manul_cond.npy')
new_data = np.load('./manul_data.npy')
new_cond_2 = np.load('./manul_cond_2.npy')
new_data_2 = np.load('./manul_data_2.npy')
new_data_total = list(new_data)
new_data_total.extend(list(new_data_2))
new_data_total = np.array(new_data_total)
new_cond_total = list(new_cond)
new_cond_total.extend(list(new_cond_2))
new_cond_total = np.array(new_cond_total)
print(new_cond_total.shape, new_data_total.shape)

# new_cond=[]
# new_data = []

# for fi in files:
#     print(fi)
#     cond = []
#     datap = []
#     with open('./new_norms/'+fi) as f:
#         reader = csv.reader(f, delimiter='\t')
#         for line in reader:
#             if len(line) == 0: #blank line
#                 continue
#             gn = line[0]
#             if line[-1] == '':
#                 line = line[:-1]
#             if gn == '' or gn not in gn2ind:
#                 if len(cond) == 0:
#                     cond = line[1:]
#                     sample_cnt = len(cond)
#                 continue
#             if len(datap) == 0:
#                 datap = np.zeros((sample_cnt, gene_len + 1), dtype='float32')
#             for i in range(sample_cnt):
#                 # print(float(cont[i+1].strip()))
#                 try:
#                     if line[i + 1].strip() == 'nan' or line[i + 1].strip() == '':
#                         datap[i][gn2ind[gn]] = 0.0
#                     else:
#                         if float(line[i + 1].strip()) == np.inf or float(line[i + 1].strip()) == -np.inf:
#                             print('----------------表达量太大')
#                             break
#                         datap[i][gn2ind[gn]] = float(line[i + 1].strip())
#                 except:
#                     print('try error')
#                     print(line)
#
#         new_data.extend(datap)
#         print(fi, np.array(datap).shape)
#         for cd in cond:
#             if not fi.startswith('GSM'):
#                 new_cond.append(fi + '---' + cd)
#             else:
#                 new_cond.append(fi.split('_')[0])

# np.save('./manul_data_2', new_data)
# np.save('./manul_cond_2', new_cond)
s2g = np.load('./manul_s2g.npy', allow_pickle=True).item()
s2g_2 = np.load('./manul_s2g_2.npy', allow_pickle=True).item()

s2g_total = {}
for sc in new_cond_total:
    if sc in s2g:
        s2g_total[sc] = s2g[sc]
    else:
        s2g_total[sc] = s2g_2[sc]
# np.save('./manul_data_total', new_data_total)
# np.save('./manul_cond_total', new_cond_total)
# np.save('./manul_s2g_total', s2g_total)
# print(new_data_total.shape, new_cond_total.shape)
assert 1==0

# # for cd in new_cond:
#     print(cd)
#     if cd.startswith('GSM'):
#         gsm = cd
#     else:
#         cnt +=1
#         gsm = input('gsm = ')
#     s2g[cd] = gsm
# np.save('manul_s2g_2', s2g)
#
assert 1==0 # alter xls data
# files = os.listdir('./norms')
# files.sort()
# for fi in files:
#     if 'xls' not in fi:
#         continue
#
#     print(fi)
#     if fi != 'GSE93232_processed_data_mRNA-seq.xlsx':
#         continue
#     reader = xlrd.open_workbook('./norms/' + fi)
#     sheet = reader.sheet_by_index(0)
#     line0 = sheet.row_values(0)
#     line1 = sheet.row_values(1)
#     print(line0)
#     print(line1)
#     dels = input('dels = ').split(' ')
#     newf = ''
#     for ri in range(sheet.nrows):
#         line = sheet.row_values(ri)
#         for i in range(len(line)):
#             if str(i) in dels:
#                 continue
#             if i == len(line) -1:
#                 newf += str(line[i])
#             else:
#                 newf += str(line[i]) + '\t'
#         newf+='\n'
#     with open('./new_norms/'+fi+'.txt', 'w') as f:
#         f.write(newf)
#
#
# assert 1==0

#process none xls file
files = os.listdir('./norms')
files.sort()
for fi in files:
    if 'xls' in fi:
        continue
    print(fi)
    if not  fi.startswith('GSM3184'):
        continue
    with open('./new_norms/' +fi, errors='ignore') as f:

        reader = csv.reader(f, delimiter='\t')
        print(fi)
        firstLine = f.readline().split('\t')
        print(firstLine)
        dels = input('dels = ').split(' ')
        newf = ''
        for c in range(len(firstLine)):
            if str(c) in dels:
                continue
            if c != len(firstLine) - 1:

                newf+= firstLine[c] + '\t'
            else:
                newf += firstLine[c]
        newf+='\n'

        for line in reader:

            for k in range(len(line)):
                if str(k) in dels:
                    continue
                if k != len(line) - 1:
                    newf += line[k] + '\t'
                else:
                    newf += line[k]
            newf += '\n'

        with open('./new_norms/' + fi, 'w') as f:
            f.write(newf)

