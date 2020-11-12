import numpy as np
import requests
import re

sample_cond = np.load('./sample_cond_new.npy')

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

def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"


print('++++++++++++++++',sample_cond[0])
test = sample_cond[0]
test_tl = []
sample_cond = sample_cond[:34]

test_old = test
sample_inf = []


for d in sample_cond:
    if not d.startswith('GSE96812'):
        continue

    print(test, '_+_+_+_+',d)

    if d.split('_')[0] == test.split('_')[0] and d != sample_cond[-1]:
        test_tl.append(d.split('---')[1])
    else:
        if d == sample_cond[-1]:
            test_tl.append(d.split('---')[1])
        gse = test.split('_')[0]

        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+gse.capitalize()
        gmslist = re.findall('GSM\d+',getHtmlText(url))
        gsmset = set(gmslist)
        gmslist = list(gsmset)
        gsm_titles =[]
        work = True
        for i in gmslist:
            # print(i)
            url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+i.capitalize()
            gsmcont = getHtmlText(url)
            ttl = re.findall('<tr valign="top"><td nowrap>Title</td>\n<td .*">([+-/A-Za-z0-9-_\s]*)</td>',gsmcont)
            # print(ttl)

            if len(ttl)>0:
                gsm_titles.append(ttl[0])
            else:
                work =False
        if not work or len(gsm_titles) <1:
            print('not work, gse = ', gse)
            test = d.copy()
            test_tl = []
            continue

        print(test_tl)
        print(gsm_titles)

        for t in test_tl:

            dists = []
            for gt in gsm_titles:
                # if gt.startswith('bpc1,2,3,4,6') and t.startswith('bpc'):
                #     gt = gt.replace('bpc1,2,3,4,6', 'bpc')
                dic = {}

                dists.append(edit_distance(t,gt))
            print(t,'---', gsm_titles[np.argmmin(dists)], '---',gmslist[np.argmmin(dists)])
            sample_inf.append(gse+'---'+t+'---'+gsm_titles[np.argmmin(dists)]+'---'+gmslist[np.argmmin(dists)])
            # gsm_titles.remove(gsm_titles[np.argmin(dists)])

        #为下一文件内的样本循环作准备

        test = d.copy()
        test_tl = []
        test_tl.append(d.split('---')[1])


# np.save('./sample_inf', sample_inf)



