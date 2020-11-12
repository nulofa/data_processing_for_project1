import numpy as np
import requests

#sam2gsm_final为添加为匹配后的，sam2gsm_final2为添加前的
s2g = np.load('./sam2gsm_final.npy', allow_pickle=True)
s2g = s2g.item()
sam_cond = np.load('./sample_cond_dels.npy')
import re

def getHtmlText(url):
    try:
        r = requests.get(url, timeout = 30)
        r.raise_for_status()
        r.encoding = r.apparent_encoding
        return r.text;
    except:
        return "产生异常"
gsm_titles = []
cnt = 0
for s in sam_cond:
    if s.startswith('GSE74488'):

        if len(gsm_titles) == 0:
            cn = s.split('---')[1]
            gse = s.split('_')[0]
            url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + gse.capitalize()
            gmslist = re.findall('GSM\d+', getHtmlText(url))
            gsmset = set(gmslist)
            gmslist = list(gsmset)
            work = True
            for i in gmslist:
                # print(i)
                url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + i.capitalize()
                gsmcont = getHtmlText(url)

                ttl = re.findall('<tr valign="top"><td nowrap>Title</td>\n<td .*">([+-/A-Za-z0-9-_\s]*)</td>',
                                 gsmcont)

                desc = re.findall('<tr valign="top"><td nowrap>Description</td>\n<td style="text-align: justify">([+-/A-Za-z0-9-_\s]*)<br></td>', gsmcont)
                print(ttl, desc)

                if len(ttl) > 0:
                    if len(desc)>0:
                        gsm_titles.append((ttl[0].replace('h','hr'),desc[0]))
                    else:
                        gsm_titles.append((ttl[0].replace('h','hr'), ''))

        for i in range(len(gsm_titles)):
            cn = s.split('---')[1]
            if 'hr' not in cn:
                if s.split('---')[1] == gsm_titles[i][1]:
                    print(s,gmslist[i])
                    cnt+=1
                    s2g[s] = gmslist[i]
            elif s.split('---')[1] == gsm_titles[i][0]:
                    print(s,gmslist[i])
                    cnt+=1
                    s2g[s] = gmslist[i]
np.save('./sam2gsm_final', s2g)
print(cnt)