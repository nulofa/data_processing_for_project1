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

sam_d = np.load('./sam_details_10_2.npy', allow_pickle=True).item()

for sc in sam_d:
    # print(sc)

    if len(sam_d[sc][1]) == 0:
        gsm = sam_d[sc][0]
        # print(sc, sam_d[sc])
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