import  numpy as np


sam2gsm = np.load('./sam2gsm.npy',allow_pickle=True)
sam2gsm =sam2gsm.item()
print(len(sam2gsm.keys()))
print(sam2gsm)
assert 1==0
sam_inf = np.load('./sample_inf.npy')

gse = ''
gse_old = ''
ttls_file = []
ttls_web = []
gsms = []
ynlist = ['y', 'n', 'y', 'y', 'n', 'y', 'y', 'y', 'y', 'n', 'y', 'n', 'y', 'y', 'n', 'n', 'y', 'y', 'y', 'n', 'n', 'n', 'n', 'n', 'y', 'y', 'n', 'n', 'n', 'n', 'y', 'y', 'n', 'n', 'y', 'n', 'n', 'y', 'n', 'n', 'y', 'n', 'y', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'y', 'n', 'y', 'n', 'n', 'n', 'n', 'y', 'n', 'n', 'n', 'n', 'n', 'n', 'y', 'n', 'n', 'n', 'y', 'y', 'n', 'n', 'y', 'y', 'y', 'y', 'n', 'y', 'y', 'n', 'n', 'y', 'y', 'n', 'y', 'n', 'y', 'y', 'n', 'n', 'y', 'n', 'n', 'n', 'n', 'n', 'y', 'n', 'n', 'y', 'n', 'n', 'n', 'n', 'n', 'n', 'y', 'y', 'n', 'n', 'n', 'y', 'y', 'n', 'n', 'y', 'y', 'n', 'n', 'y', 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'y', 'n', 'n', 'n', 'n']
i = 0
j = 0
i4yn = 0
for si in sam_inf:
    # print(si)
    if gse == '':
        gse = si.split('---')[0]
        gse_old = gse

        ttls_file.append(gse+'---'+si.split('---')[1])
        ttls_web.append(gse+'---'+si.split('---')[2])
        gsms.append(si.split('---')[3])
        i +=1
    else:
        gse = si.split('---')[0]
        if gse == gse_old:

            ttls_file.append(gse+'---'+si.split('---')[1])
            ttls_web.append(gse+'---'+si.split('---')[2])
            gsms.append(si.split('---')[3])
            i+=1
        else:

            gse_old =gse
            for k in range(j,i):
                print(ttls_file[k].split('---')[1], '------', ttls_web[k].split('---')[1], '--------',gsms[k])
            # flag = input('y or n: ')
            flag = ynlist[i4yn]
            i4yn+=1

            if flag == 'n':
                i=j
                ttls_file = ttls_file[:j]
                ttls_web = ttls_web[:j]
                gsms = gsms[:j]
                # print('n:------')
                # for k in range(len(ttls_file)):
                #     print(ttls_file[k], ttls_file[k])
                # print('n:------end')
            else:
                j = i
                # print('y:------')
                # for k in range(len(ttls_file)):
                #     print(ttls_file[k], ttls_file[k])
                # print('y:------end')
            ttls_file.append(gse+'---'+si.split('---')[1])
            ttls_web.append(gse+'---'+si.split('---')[2])
            gsms.append(si.split('---')[3])
            i+=1
sam2gsm ={}
for k in range(len(ttls_file)):
    sam2gsm[ttls_file[k]] = gsms[k]
    # print(sam2gsm)
np.save('./sam2gsm', sam2gsm)