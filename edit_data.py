import csv
import os

import numpy as np
path = './data0/'
path2 = './data1/'
file = 'GSE43703_processed_data_repeat2_samples.txt'

data = []
dataset = []
# files = os.listdir(path)
# for file in files:

data = []
with open(path+file, encoding='utf-8') as f:
    if file.endswith('csv'):
        f_csv = csv.reader(f)
        cont = ''
        for line in f_csv:
            ds = '\t'.join([d for d in line[1:]])
            cont = line[0].split('.')[0] + '\t' + ds
            data.append(cont)
    else:

        for line in f:

            cont = line.split('\t')

            # cont = line.split('\t')
            # cont[2] = cont[2].replace('"','')
            # cont[0] =cont[0].replace('"','')
            new_cont = '\t'
            if cont[0] != '':
                new_cont = cont[0].split('.')[0]+'\t'


            for c in cont[16:17]:

                new_cont+=c.rstrip()+'\t'
                # print(new_cont)
            new_cont.rstrip()
            data.append(new_cont)


# np.save('iso_samples', samples)
# np.save('iso_vals', vals)
# np.save('iso_genes', genes)
# dataset = np.array(dataset).T
#
# with open(path+file2, 'w') as f:
#     for data in dataset:
#         dlen = len(data)
#         i = 0
#         for d in data:
#             i+=1
#             if i == dlen:
#                 f.write(d+'\n')
#             else:
#                 f.write(d+'\t')


    with open(path2+file.replace('.csv', '.txt'), 'w') as f:
        for d in data:
            f.write(d+'\n')
    # #         # print(d)
    # #         # assert 1==0

