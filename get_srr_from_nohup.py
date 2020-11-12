import re

srrlist = []
with open('./nohup') as f:
    for line in f:
        srr = re.findall('failed to download (SRR.*)', line)
        if len(srr) > 0:
            srrlist.append(srr[0])
print(srrlist)