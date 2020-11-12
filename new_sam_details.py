import numpy as np
import requests
import re

def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"


# ttls = ['cufflinks_ath_leaf_rep1','cufflinks_ath_leaf_rep2','cufflinks_ath_leaf_rep3','cufflinks_ath_buds_rep1','cufflinks_ath_buds_rep2','cufflinks_ath_buds_rep3']
# gmslist = ['GSM2046083n','GSM2834608','GSM2834610','GSM2834612','GSM3508146GSM3508147']
sam2gsm = np.load('./sam2gsm_final.npy',allow_pickle=True)
sam2gsm = sam2gsm.item()
sam_details = np.load('./sam_details_dels2.npy',allow_pickle=True)
sam_details=sam_details.item()
i=0
gts=[]
dels =[]
for gt in sam2gsm:
    # if len(sam_details[gt][1])>0:
    #     continue

    gsm = sam2gsm[gt]
    gsm = gsm.strip()
    print(gsm)
    # if gsm == 'GSM2046083n':
    #     dels.append(gt)
    # if gsm == 'GSM3508146GSM3508147':
    #     gsm = 'GSM3508147'
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+gsm
    gsmcont = getHtmlText(url)
    # tissue = re.findall('tissue: ([^:]*)<br>', gsmcont)
    title = re.findall('<tr valign="top"><td nowrap>Title</td>\n<td style="text-align: justify">([^\n]*)</td>\n</tr>',gsmcont)

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
    characteristics = re.findall('<tr valign="top"><td nowrap>Characteristics</td>\n<td style="text-align: justify">(.*)<br></td>', gsmcont)
    # sam_details[gse+ttl] = (gsm, tissue, genotype, treatment)

    # print((gsm, tissue, genotype, treatment))
    sam_details[gt] = (gsm, title, characteristics)
    print((gsm, title, characteristics), gt)

    i+=1

for d in dels:
    sam_details.pop(d)
# np.save('./sam_details_dels2', sam_details)

