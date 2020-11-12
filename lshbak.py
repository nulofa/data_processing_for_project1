import numpy as np
import math
import time
import os
import csv

import matplotlib.pyplot as plt

def read_data0(mode):
    if mode == 'a':
        dataset = np.load('./dataset.npy')
        sample_id = np.load('./sample_id.npy')
        sample_cond = np.load('./sample_cond.npy',allow_pickle=True)
        data_num = len(dataset)
        dim = len(dataset[0])
        return dataset,data_num, dim, sample_id, sample_cond

    gene_name = np.load('./gene_name.npy')
    gn2ind = {}
    gene_len = len(gene_name)
    for i in range(gene_len):
        gn2ind[gene_name[i]] = i

    path = './data/'
    files = os.listdir(path)
    print(files)
    file_cnt = 0
    sample_id = []
    sample_cond = {}
    cnt = 0
    for file in files:
        file_cnt+=1
        print(file)

        if file_cnt ==10:
            break
        with open(path+file) as f:
            flag = False
            ind = 0
            log_flag = False
            for line in f:
                if line.startswith('!dataset_value_type'):
                    if line.split('=')[1].strip()[:3] == 'log':
                        print('log ratio ',file)
                        log_flag = True
                if line.startswith('#GSM'):
                    sampleid = line.split('=')[0].rstrip() #包含符号#
                    sample_id.append(sampleid[1:])
                    condition = line.split(':')
                    sample_cond[sampleid[1:]] = condition[1]+':'+condition[2].rstrip()

                if line.startswith('!dataset_sample_count'):
                    sample_count = int(line.split('=')[1])

                if line.startswith('ID_REF') :
                    flag = True
                    data = np.zeros((gene_len, sample_count))
                elif flag:  #表达数据开始
                    if line.startswith('!dataset_table_end'):
                        break
                    gn = line.split('\t')[1]
                    string = line.split('\t')[2:2+sample_count]
                    if gn in gn2ind:
                        ind = gn2ind[gn]
                    else:
                        cnt +=1
                        if gn.strip() in gene_name:
                            print('fffffffffffffffffffffff')
                        if gn == '--control' or gn =='--Control':
                            break
                        else:
                            continue

                    for i in range(sample_count):
                        if string[i].rstrip() != 'null':
                            data[ind][i] = float(string[i])
                        else:
                            data[ind][i] = 0

                    zero_ind = data[ind] == 0
                    if not log_flag:
                        avg = np.average(data[ind])
                        if np.sum(data[ind] <0) > 0 or avg ==0:    #has negetive value
                            # print(avg+eps)
                            up_ind = data[ind]>avg
                            data[ind] = [-1 for i in range(sample_count)]
                            data[ind][up_ind]=1
                            data[ind][zero_ind]=0
                        else:
                            data[ind][zero_ind] = avg
                            data[ind] = np.log2(data[ind] / avg)

        dataset = data
        if file_cnt ==1:
            data_set = dataset
        else:
            data_set = np.hstack((data_set,dataset))
    # print(length)
    data_set = data_set.T
    print(len(data_set), len(data_set.T))
    print('cnt = ',cnt)
    print(gn2ind)
    return data_set,len(data_set),gene_len, sample_cond

def read_data(mode):
    if mode == 'a':
        data_set = np.load('./dataset.npy')
        sample_cond = np.load('./sample_cond.npy')
        sample_dup = np.load('./sample_dup.npy', allow_pickle=True)
        return data_set, data_set.shape[0], data_set.shape[1],sample_cond, sample_dup

    #data0 是exp和fpkm的搜索结果
    path = './data0/'
    files = os.listdir(path)
    txt_files = []
    csv_files = []
    sample_dup = {}
    for fs in files:
        if fs.endswith('.txt') :
            txt_files.append(fs)
        if fs.endswith('.csv'):
            csv_files.append(fs)

    print(files)
    sample_cond = []
    file_cnt = 0
    gene_name = np.load('./gene_name0.npy')
    gn_dup = np.load('./gn_dup.npy',allow_pickle=True)

    gn_dup = gn_dup.item() #np 读取保存字典必须
    gn2ind = {}
    gene_len = len(gene_name)
    print(gene_len)
    cnt = 0
    for i in range(gene_len):
       gn2ind[gene_name[i]] = i
       if gene_name[i] in gn_dup:
           for gn in gn_dup[gene_name[i]]:
               gn2ind[gn] = i

    data_set = []
    # txt_files=txt_files[:20]
    for file in txt_files:
        file_cnt+=1
        # print(file)
        if 'counts' in file.lower():
            print('count file = ', file)
            continue
        cond = []
        with open(path+file,errors='ignore') as f:
            init_datap = False
            exp_file = False
            break_flag = False
            enter_once = False
            log_flag = False
            if 'log' in file.lower():
                log_flag = True
            for line in f:
                if break_flag:
                    break
                cont = line.split('\t')
                gn = cont[0].strip()
                if file == 'GSE107018_fpkm.txt':
                    cont = cont[1:]
                # sample_cnt = len(cont) - 1
                # if len(cond)>0 and len(cond) != sample_cnt: #cond 已经初始化，但与cnt不同
                #     print(cond, sample_cnt)
                #     print('file ---------------', file)
                #     break
                if gn == '' or gn not in gn2ind:
                    if len(cond) == 0:
                        cond.extend(cont[1:])
                        sample_cnt = len(cont) - 1
                        # print(cond)
                    continue

                enter_once = True
                # print(file, sample_cnt, line.split('\t')[1:])
                if not init_datap:
                    init_datap = True
                    datap = np.zeros((sample_cnt,gene_len))

                for i in range(sample_cnt):
                   # print(float(cont[i+1].strip()))
                   try:
                       datap[i][gn2ind[gn]] = float(cont[i+1].strip())
                   except:
                       # print('except, ', cont[i+1].strip())
                       break_flag = True
                       break

                zero_ind = datap[:,gn2ind[gn]] == 0

                # if not log_flag:
                #
                #     # assert 1==0
                #     avg = np.average(datap[:,gn2ind[gn]])
                #
                #     if np.sum(datap[:,gn2ind[gn]] < 0) > 0 or avg == 0:  # has negetive value
                #        # print(avg+eps)
                #        # pass
                #        up_ind = datap[:,gn2ind[gn]] > avg
                #        datap[:,gn2ind[gn]] = [-1 for i in range(sample_cnt)]
                #        datap[up_ind,gn2ind[gn]] = 1
                #        datap[zero_ind,gn2ind[gn]] = 0
                #     else:
                #        # if (np.sum(datap[:, gn2ind[gn]]) - np.max(datap[:, gn2ind[gn]]))*5 < np.max(datap[:, gn2ind[gn]]):
                #        #      print('before: ', datap[:, gn2ind[gn]], file)
                #        # print('avg = ', avg)
                #        # print(datap[zero_ind,gn2ind[gn]] )
                #        datap[zero_ind,gn2ind[gn]] = avg
                #        datap[:,gn2ind[gn]] = np.log2(datap[:,gn2ind[gn]] / avg)


            else: #for..else, 上面for正常才执行else
                if enter_once:
                    exp_file = True
        if not exp_file:
            print('not standard exp : ', file)
        else:
            data_set.extend(datap)

            sample_id = file.split('.')[0]
            last_cd = ''
            for cd in cond:
                if cd[:-1] == last_cd[:-1]:
                    if cd[:-1] not in sample_dup:
                        sample_dup[cd[:-1]] = [last_cd]
                    else:
                        if last_cd not in sample_dup[cd[:-1]]:
                            sample_dup[cd[:-1]].append(last_cd)
                    if cd not in sample_dup[cd[:-1]]:
                        sample_dup[cd[:-1]].append(cd)
                sample_cond.append(sample_id+'---'+cd.rstrip())
                last_cd = cd
                if cd == '':
                    print('empt = ',file)
            print(file, sample_cnt, len(sample_cond))


    # data_set = []
    log_flag = False

    for file in csv_files:
        if 'log' in file.lower():
            log_flag = True
        if 'counts' in file.lower():
            continue
        with open(path+file) as f:
            f_csv = csv.reader(f)
            sample_cnt = -1
            datap = []
            break_flag = False
            enter_once = False
            for line in f_csv:
                if sample_cnt == -1:
                    sample_cnt = len(line)-1
                    cond= line[1:]
                elif len(datap) == 0:
                    datap = np.zeros((sample_cnt, gene_len))
                else:
                    gn = line[0].strip()
                    if gn not in gn2ind:
                        # print(gn, 'not in gn2ind')
                        continue
                    enter_once = True
                    for i in range(sample_cnt):
                        try:
                            datap[i][gn2ind[gn]] = float(line[i + 1].strip())
                        except:
                            print('csv format error : ',file)
                            break_flag = True
                            break
                    if break_flag == True:
                        break

                    zero_ind = datap[:, gn2ind[gn]] == 0

                    # if not log_flag:
                    #
                    #     # assert 1==0
                    #     avg = np.average(datap[:, gn2ind[gn]])
                    #
                    #     if np.sum(datap[:, gn2ind[gn]] < 0) > 0 or avg == 0:  # has negetive value
                    #         # print(avg+eps)
                    #         # pass
                    #         up_ind = datap[:, gn2ind[gn]] > avg
                    #         datap[:, gn2ind[gn]] = [-1 for i in range(sample_cnt)]
                    #         datap[up_ind, gn2ind[gn]] = 1
                    #         datap[zero_ind, gn2ind[gn]] = 0
                    #     else:
                    #         # if (np.sum(datap[:, gn2ind[gn]]) - np.max(datap[:, gn2ind[gn]]))*5 < np.max(datap[:, gn2ind[gn]]):
                    #         #      print('before: ', datap[:, gn2ind[gn]], file)
                    #         # print('avg = ', avg)
                    #         # print(datap[zero_ind,gn2ind[gn]] )
                    #         datap[zero_ind, gn2ind[gn]] = avg
                    #         datap[:, gn2ind[gn]] = np.log2(datap[:, gn2ind[gn]] / avg)

            else:#for...else
                if enter_once:
                    data_set.extend(datap)
                    for cd in cond:

                        if cd[:-1] == last_cd[:-1]:
                            if cd[:-1] not in sample_dup:
                                sample_dup[cd[:-1]] = [last_cd]
                            else:
                                if last_cd not in sample_dup[cd[:-1]]:
                                    sample_dup[cd[:-1]].append(last_cd)
                            if cd not in sample_dup[cd[:-1]]:
                                sample_dup[cd[:-1]].append(cd)
                        last_cd = cd
                        if cd == '':
                            print('empt = ', file)
                        sample_cond.append(sample_id + '---' + cd.rstrip())
                    print(file, sample_cnt, len(sample_cond))
                else:
                    print('not exp csv：', file)

    data_set = np.array(data_set)
    print(data_set.shape)
    # print((sample_dup))
    # assert 1==0
    return data_set, data_set.shape[0], data_set.shape[1],sample_cond, sample_dup

def signature_bit(data, planes):
    """
    LSH signature generation using random projection
    Returns the signature bits for two data points.
    The signature bits of the two points are different
     only for the plane that divides the two points.
     """
    sig = 0
    for p in planes:
        sig <<= 1 #左移一位，相当于乘2
        if np.dot(data, p) >= 0:
            sig |= 1 #点在平面下方，该位为1
    return sig

def length(v):
    """returns the length of a vector"""
    return math.sqrt(np.dot(v, v))

def l2(v1,v2):
    assert len(v1) == len(v2)
    v1 =np.array(v1)
    v2 =np.array(v2)
    return np.sum((v1-v2)**2)**0.5

def pearson(vector1, vector2):
    n = len(vector1)
    #simple sums
    # sum1 = sum(float(vector1[i]) for i in range(n))
    # sum2 = sum(float(vector2[i]) for i in range(n))
    # #sum up the squares
    # sum1_pow = sum([pow(v, 2.0) for v in vector1])
    # sum2_pow = sum([pow(v, 2.0) for v in vector2])
    # #sum up the products
    # p_sum = sum([vector1[i]*vector2[i] for i in range(n)])

    sum1=0;sum2=0;sum1_pow=0;sum2_pow=0;p_sum=0
    for i in range(n):
        sum1 += vector1[i]
        sum2 += vector2[i]
        sum1_pow += vector1[i]**2.0
        sum2_pow += vector2[i] ** 2.0
        p_sum += vector1[i]*vector2[i]
    #分子num，分母den
    num = p_sum - (sum1*sum2/n)
    den = math.sqrt((sum1_pow-pow(sum1, 2)/n)*(sum2_pow-pow(sum2, 2)/n))
    if den == 0:
        return 0.0
    return num/den


def cos_sim(v1, v2):

    assert len(v1) == len(v2)
    # return pearson(v1,v2)
    return np.dot(v1,v2)/ length(v1)/length(v2)

def linear_search(data, q):

    total_sim = []
    max_cos= cos_sim(data[0], q)
    res = 0
    data_num = len(data)
    for p in range(data_num):
        cosine = cos_sim(data[p], q)
        total_sim.append(cosine)
        if cosine > max_cos:
            max_cos = cosine
            res = p
        # print('cos = ',cosine)
    return res, max_cos,total_sim

def val_to_rank1(data):
    dataset = data.copy()
    dataset = np.array(dataset)
    data_len = len(dataset)
    dim = len(dataset[0])

    for i in range(data_len):
        ind = np.argsort(dataset[i])
        ind2 = dataset[i] == (np.zeros((1,dim)))[0]
        dataset[i][ind] = range(1,dim+1)
        dataset[i][ind2] = 0
    return np.array(dataset)

def lsh_tables(data, bits, dim,iter):
    tables = []
    planes = []
    data_num = len(data)
    for i in range(iter):
        ref_planes = np.random.randn(bits, dim)
        planes.append(ref_planes)
        table = {}

        for r in range(data_num):
            # signature bits for data points
            sig = signature_bit(data[r], ref_planes)

            if sig not in table:
                table[sig] = [r]
            else:
                table[sig].append(r)

        tables.append(table)
    return tables, planes

def roc_curve(dataset):
    xs = []
    ys = []
    #val_sample_ind 没有用
    dataset = np.array(dataset)
    dim = dataset.shape[1]
    data_num = dataset.shape[0]
    shift = 0
    triples_num = data_num//3
    for i in range(triples_num):

        for ind in range(3):
            query_q = np.zeros(dim)
            labels = np.zeros((data_num, 1))
            query_q += dataset[ind+shift*3]
            # q_ind = np.argsort(query_q)
            # query_q[q_ind] = range(1, dim + 1)
            # zero_ind = query_q ==(np.zeros((1,dim)))[0]
            # query_q[zero_ind] = 0
            if ind == 0:
                labels[shift*3+ind:shift*3+ind+3] =1
            if ind ==1:
                labels[shift * 3 + ind] = 1
                labels[shift * 3 + ind-1] = 1
                labels[shift * 3 + ind+1] = 1
            if ind ==2:
                labels[shift * 3 + ind] = 1
                labels[shift * 3 + ind-1] = 1
                labels[shift * 3 + ind-2] = 1
            # print('ind = ', ind + shift * 3, 'label = ', labels)

            res1,res2,sim = linear_search(dataset,query_q)

            arg_ind = np.argsort(sim)
            x=[0]
            y=[0]
            arg_ind =arg_ind[::-1]

            for ind in arg_ind:
                i = sim>= sim[ind]
                j = sim <sim[ind]
                tp = sum(labels[i])
                fn = sum(labels[j])
                tn = sum(labels[j]==0)
                fp = sum(labels[i]==0)
                x.append(fp/(fp+tn))
                y.append(tp/(tp+fn))
            if len(xs)==0:
                xs = np.array(x)
            else:
                xs +=x
            if len(ys) == 0:
                ys = np.array(y)
            else:
                ys += y

        shift += 1 #下一个3重复组

    xs = xs/triples_num/3
    ys = ys/triples_num/3

    # print(xs[-1],ys[-1])
    return xs, ys

if __name__ == '__main__':
    # data = np.random.rand(3,4)+1
    # data[data<1.5] = 0
    # avg = np.average(data,axis=0)
    # avg[avg ==0] =0.01
    # print(data)
    # for d in data:
    #     ind = d==0
    #     d[ind] = avg[ind]
    # print(np.log(data/avg))
    # assert 1 == 0
    mode= 'a'
    dataset,data_num, dim, sample_cond, sample_dup= read_data(mode)
    # avg = np.median(dataset, axis=0)
    # avg[avg <= 0] = 0.01
    # for d in dataset:
    #     ind = d == 0
    #     if np.sum(d<0) >0:
    #         up_ind = d>avg
    #         d[:] = -1
    #         d[up_ind]=1
    #         d[ind] =0
    #     else:
    #         d[ind] = avg[ind]
    #         d[:] = np.log(d/avg)
    # dataset = val_to_rank1(dataset)
    bits =10
    table_num = 5
    if mode == 'a':
        sample_dup = sample_dup.item()
    print(dataset.shape)
    skip = False
    val_sample = []

    for sd in sample_dup.items():
        if len(sd[1]) <4 and len(sd[1]) >2:
            if sd[0] == 'FPKM.3S':
                skip = False
            if not skip:
                if sd[1][0] != '': #有一个键为空的元素
                    tmp = []
                    tmp.extend(sd[1])
                    val_sample.append(tmp)
                if sd[0] == 'wol_cut0h_': #跳过wol_cut0h_之后，FPKM.3S之前的项
                    skip = True
    print(val_sample)

    sample2ind = {}

    for i in range(len(sample_cond)):
        sample = sample_cond[i].split('---')[-1]
        sample2ind[sample] = i

    val_sample_ind = []

    for s in val_sample:
        sample_ind = []
        for s0 in s:
            sample_ind.append(sample2ind[s0])
        val_sample_ind.append(sample_ind)

    print(val_sample_ind)
    # np.save('./sample_cond', sample_cond)
    # np.save('./sample_dup', sample_dup)
    # # np.save('./dataset',dataset)
    # assert 1==0


    # assert 1==0
    # dataset = val_to_rank1(dataset)
    # val_sample_ind = val_sample_ind[:10]
    xs=[]
    ys=[]
    count =0
    val_inds = []
    for i in val_sample_ind:
        for j in i:
            val_inds.append(j)

    val_dataset = []
    for vi in val_inds:
        val_dataset.append(dataset[vi])
    x, y = roc_curve(val_dataset)

    plt.plot(x,y, '--')
    plt.show()
    assert 1==0
    # ind4q = np.random.randint(0, data_num)

    query_q = np.zeros(dim)
    true_labels = np.zeros((data_num,1))
    for ind in val_sample_ind[5]:
        query_q += dataset[ind]
        true_labels[ind] = 1
    query_q = query_q/len(val_sample_ind[5])
    # dataset = np.delete(dataset, ind4q, axis=0)
    print('query_q:', query_q.shape,sample_cond[val_sample_ind[5][0]])
    # sample_cond = np.delete(sample_cond, ind4q)
    tables, planes = lsh_tables(dataset, bits, dim, table_num)
    # print(tables)

    # for i in range(data_num):
    #     num_of_zeros = np.sum(dataset[i][:]==np.zeros((1,dim)))
    #     print(dim - num_of_zeros)





    query_num = 1
    cnt = 0
    dif = 0
    time1 = 0
    time2 = 0
    for qn in range(query_num):
        # ind1 = np.random.randint(0,data_num)
        # ind2 = np.random.randint(0, data_num)
        # # print(ind1,ind2)
        #
        # query_q =dataset[ind1].copy()
        # # print('query_q:', ind1 ,sample_id[ind1], sample_cond[sample_id[ind1]])
        #
        # index = [i%2==0 for i in range(dim)]
        # query_q[index] = dataset[ind2][index]


        # query_q = val_to_rank([query_q])[0]
        s1 = time.clock()
        total = []
        for i in range(table_num):

            sig4q = signature_bit(query_q, planes[i])

            bins = []
            bit4q = bin(sig4q)
            bit4q = bit4q[2:]
            bins.append(int(bit4q, 2))
            # print('bit4q= ',bit4q)
            for b in range(len(bit4q)):
                if bit4q[b] == '0':
                    # print(bit4q[:b]+'1'+bit4q[b+1:])
                    bins.append(int(bit4q[:b] + '1' + bit4q[b + 1:], 2))
                else:
                    bins.append(int(bit4q[:b] + '0' + bit4q[b + 1:], 2))


            for b in bins:
                if b not in tables[i]:
                    continue
                ind = tables[i][b][0]
                # print('len = ',len(tables[i][b]))
                max_cos = cos_sim(query_q, dataset[ind])
                res = ind

                for p in tables[i][b]:
                    cosine = cos_sim(dataset[p], query_q)
                    if cosine > max_cos:
                        max_cos = cosine
                        res = p
                # print(max_cos, res[:10])
                total.append((max_cos, res))
        if len(total) == 0:
            print('not found-----------------')
            continue
        res1 = max(total)
        # res1 = np.argsort(np.array(total)[:,0])
        e1 = time.clock()

        id = int(res1[1])

        condition = sample_cond[id]
        print(res1, e1-s1, 'cond : '+condition)

        s2 = time.clock()
        res2, cos2,total_sim = linear_search(dataset, query_q)
        e2 = time.clock()

        id = int(res2)

        condition = sample_cond[id]
        print(cos2, res2, e2 - s2, 'cond : '+condition)
        print(dataset[res2][:10])
        # print((e1-s1)/(e2-s2))
        if res1[1] == res2:
            cnt += 1

        time1 +=e1-s1
        time2 +=e2-s2

    print(cnt / query_num)
    print(time1/time2)


