import re
import numpy as np
import requests

sam_c = np.load('./sample_cond_dels2.npy')

gselist = []

def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;top

    except:
        return "产生异常"

for sc in sam_c:
    gse = sc.split('_')[0]
    if gse not in gselist:
        gselist.append(gse)

with open('./gds_result.txt', encoding='utf-8') as f:
    gsetotal = re.findall('Accession: (GSE[0-9]*)', f.read())
cnt = 0
for g in gsetotal:
    flag = False
    if g not in gselist:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + g.capitalize()
        # gsmlist = re.findall('GSM\d+', getHtmlText(url))
        # gsmset = set(gsmlist)
        # gsmlist = list(gsmset)
        # print(len(gsmlist))
        # cnt += len(gsmlist)
        bioproject = re.findall('">(SR.*)</a></td>', getHtmlText(url))
        if len(bioproject) ==0:
            print(url)
            print('跳过————————————————————————————————')
            continue
        biop_url = 'https://www.ncbi.nlm.nih.gov/Traces/study/?acc=' + bioproject[0]
        print(biop_url)

