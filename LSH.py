import heapq
import pickle
from scipy.stats.stats import pearsonr,spearmanr
import nanopq
import numpy as np
import math
import time
import os
import csv
from sklearn import metrics
import matplotlib.pyplot as plt

def read_data0(mode):
    if mode == 'a':
        dataset = np.load('./dataset_float16.npy')
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
        x1 = np.load('./x1.npy')  #  1:euclid 2:cos 3:pearson 4:spearman 5:lsh 6:pq
        y1 = np.load('./y1.npy')
        x2 = np.load('./x2.npy')
        y2 = np.load('./y2.npy')
        x3 = np.load('./x3.npy')
        y3 = np.load('./y3.npy')
        x4 = np.load('./x4.npy')
        y4 = np.load('./y4.npy')
        x5 = np.load('./x5.npy')
        y5 = np.load('./y5.npy')
        x6 = np.load('./x6.npy')
        y6 = np.load('./y6.npy')
        a1 = metrics.auc(x1, y1)
        a2 = metrics.auc(x2, y2)
        a3 = metrics.auc(x3, y3)
        a4 = metrics.auc(x4, y4)
        a5 = metrics.auc(x5, y5)
        a6 = metrics.auc(x6, y6)
                                # x1 = np.insert(x1,0, values=0,axis=0)
                                # y1 = np.insert(y1, 0, values=0, axis=0)
                                # x2 = np.insert(x2, 0, values=0, axis=0)
                                # y2 = np.insert(y2, 0, values=0, axis=0)
                                # x3 = np.insert(x3, 0, values=0, axis=0)
                                # y3 = np.insert(y3, 0, values=0, axis=0)
                                # x4 = np.insert(x4, 0, values=0, axis=0)
                                # y4 = np.insert(y4, 0, values=0, axis=0)   #5,6 不需要
        # plt.plot(x1,y1, '-o')
        # plt.plot(x2, y2, '-o')
        # plt.plot(x3, y3, '-o')
        plt.plot(x4, y4, '-o')
        plt.plot(x5, y5, '-o')
        plt.plot(x6, y6, '-o')

        # plt.legend(['euclidean: ' + str(round(a1, 3)), 'cosine: ' +str(round(a2,3)), 'pearson: '+str(round(a3,3))
        #            ,'spearman: '+str(round(a4,3))])
        # plt.savefig('./roc0_10_21.jpg')
        plt.legend(['spearman: '+ str(round(a4, 3)), 'pq: ' + str(round(a5, 3)), 'lsh: ' + str(round(a6, 3))])
        plt.savefig('./roc1_10_21.jpg')
        plt.show()
        assert 1==0
        data_set = np.load('./dataset_new_dels2.npy')
        sample_cond = np.load('./sample_cond_dels2.npy')
        sample_dup = np.load('./sample_dup_new.npy', allow_pickle=True)
        # ds = []
        # for d in data_set:
        #     ds.append([d])
        # data_set = np.array(ds)
        return data_set, data_set.shape[0], data_set.shape[1],sample_cond, sample_dup

    #data0 是exp和fpkm,rpkm,tpm,rpm的搜索结果
    #data1,2是测试目录
    path = './data/'
    files = os.listdir(path)
    txt_files = []
    csv_files = []
    sample_dup = {}
    for fs in files:
        if fs.endswith('.tsv') or fs.endswith('txt'):
            txt_files.append(fs)
        if fs.endswith('.csv') :
            csv_files.append(fs)

    print(files)
    sample_cond = []
    file_cnt = 0
    gene_name = np.load('./gene_name2.npy')
    gn_dup = np.load('./gn_dup2.npy',allow_pickle=True)

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
        if 'count' in file.lower() or 'reads' in file.lower() or 'raw' in file.lower():
            print('count file = ', file)
            continue
        # if 'library' in file.lower():
        #     print('library count file = ', file)
        #     continue
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
                line = line.rstrip()
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
                        # print(sample_cnt)
                    continue

                enter_once = True
                # print(file, sample_cnt, line.split('\t')[1:])
                if not init_datap:
                    init_datap = True
                    datap = np.zeros((sample_cnt,gene_len),dtype='float32')

                for i in range(sample_cnt):
                   # print(float(cont[i+1].strip()))
                   try:
                       if cont[i + 1].strip() == 'nan':
                           datap[i][gn2ind[gn]] = 0.0
                       else:
                           if float(cont[i+1].strip()) == np.inf or float(cont[i+1].strip()) == -np.inf:
                               print('----------------表达量太大', file)
                               break_flag = True
                               break
                           datap[i][gn2ind[gn]] = float(cont[i + 1].strip())


                   except:
                       print('except, --------', cont)
                       break_flag = True
                       break

                # zero_ind = datap[:,gn2ind[gn]] == 0

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
            # if '副本' not in file:
            #     os.remove(path+file)
        else:
            data_set.extend(datap)

            sample_id = os.path.splitext(file)[0]
            if '副本2' in sample_id:
                sample_id = sample_id.replace('副本2','')
            if '副本' in sample_id:
                sample_id = sample_id.replace('副本','')
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
                if line[-1] == '':
                    line = line[:-1]
                if sample_cnt == -1:
                    sample_cnt = len(line)-1
                    cond= line[1:]
                elif len(datap) == 0:
                    datap = np.zeros((sample_cnt, gene_len),dtype='float32')
                else:
                    gn = line[0].split(',')[0].strip()

                    if gn not in gn2ind:
                        # print(gn, 'not in gn2ind')
                        continue
                    enter_once = True
                    for i in range(sample_cnt):
                        try:
                            if float(line[i + 1].strip()) == np.inf:
                                print('----------------表达量太大', file)
                                exit('end')
                                break_flag = True
                                break
                            datap[i][gn2ind[gn]] = float(line[i + 1].strip())
                        except:
                            print('csv format error : ',file, line[i + 1].strip())
                            print(line)
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
                    sample_id = os.path.splitext(file)[0]
                    if '副本' in sample_id:
                        sample_id = sample_id.replace('副本', '')
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
    # print(type(v))
    v = v.astype('float64')
    return math.sqrt(np.dot(v, v))

def l2(v1,v2):
    # return hamming_dist(v1, v2)

    assert len(v1) == len(v2)
    v1 =np.array(v1).astype('float32')
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

    # v1 = np.array(v1).reshape(1,-1)
    # v2 = np.array(v2).reshape(1,-1)
    # print(length(v1), length(v2))
    assert len(v1) == len(v2)
    # return pearson(v1,v2)
    # v1 = v1.astype('float32')
    # v2 = v2.astype('float32')
    # print('dot', np.dot(v1,v2),length(v1),np.dot(v1,v1))
    return np.dot(v1,v2)/ length(v1)/length(v2)

def count1s(m):
    cnt = 0
    while(m):

        m = m & (m-1)
        cnt += 1
    return cnt

def hamming_dist(v1, v2):
    # print(v1)
    # print(v2)

    assert len(v1) == len(v2)
    l = 162
    # a = v1[0]
    # b = v2[0]
    # print(a,b,v2[0].dtype)
    # assert 1==0
    # print(count1s(v1[0]^v2[0]))
    sum1s = count1s(v1[0]^v2[0])    #no. of not identical bit
    return sum1s/l

def linear_search(data, q, k, ri):

    total_sim = []
    max_cos= -1
    # max_cos = hamming_dist(data[0], q)

    res = 0
    data_num = len(data)
    topk_sim = [(0.0, -1) for i in range(k)]
    heapq.heapify(topk_sim)
    for p in range(data_num):
        if ri ==0:
            cosine = pearsonr(data[p], q)[0]
            # print(cosine)
        # print(cosine, ' cos')
        elif ri==1:
            cosine = l2(data[p], q)


        elif ri==2:
            cosine = spearmanr(data[p], q)[0]

        else:
            cosine = cos_sim(data[p], q)
        # cosine = hamming_dist(data[p], q)   #for lsh
        total_sim.append(cosine)
        if ri != 1:
            if cosine > topk_sim[0][0] and np.sum([p == i[1] for i in topk_sim]) == 0:
                # print('----before: ',topk_sim, cosine)
                heapq.heappop(topk_sim)
                heapq.heappush(topk_sim, (cosine, p))
                # print('-----after: ', topk_sim)
            if cosine > max_cos:
                max_cos = cosine
                res = p
        else:
            pass
        # print('cos = ',cosine)

    return res, max_cos,total_sim, topk_sim

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
        if i == 0:
            data_lsh = []
        for r in range(data_num):
            # signature bits for data points
            sig = signature_bit(data[r], ref_planes)
            if i ==0:
                data_lsh.append([sig])
            if sig not in table:
                table[sig] = [r]
            else:
                table[sig].append(r)
        if i ==0:
            np.save('./data_lsh_104', np.array(data_lsh))
        tables.append(table)
    return tables, planes


def linear_search2(dataset, query_q,pq):
    dists = pq.dtable(query=query_q).adist(codes=dataset)
    # print(dists, len(dists))
    return dists

def roc_curve2(dataset, ri,pq, k_dup):
    xs = []
    ys = []
    #val_sample_ind 没有用
    # dataset = np.array(dataset, dtype='int64')

    dataset = np.array(dataset)


    data_num = dataset.shape[0]
    shift = 0
    triples_num = data_num//k_dup

    for i in range(triples_num):

        for ind in range(k_dup):
            # query_q = np.zeros((1,27))
            # query_q = np.zeros(dim)


            # if ind != 0:
            #     continue
            query_q = dataset[ind+shift*k_dup].copy()
            # print(query_q)
            labels = np.zeros((data_num, 1))

            labels[shift * k_dup:shift * k_dup + k_dup] = 1
            # print(ind+shift*k_dup,labels[:21])
            if pq != None:
                print('pq:--')
                #asymmetric
                # dataset2 = pq.encode(dataset)
                # query_q = dataset[ind+shift*k_dup]

                #symmetric
                dataset2 = dataset.reshape(2*k_dup,-1)
                query_q = pq.decode(dataset[ind+shift*k_dup])[0]
                sim = linear_search2(dataset2, query_q, pq)
            else:
                res1,res2,sim ,topk_sim= linear_search(dataset,query_q,1,ri)
            sim = np.array(sim)
            arg_ind = np.argsort(sim)


            x=[0]
            y=[0]
            if ri != 1:
                arg_ind =arg_ind[::-1]

            top_n = 20
            sim = sim[arg_ind[:top_n]]
            # print(sim)
            labels = labels[arg_ind[:top_n]]
            # print(labels)
            # fpr, tpr, thresholds = metrics.roc_curve(labels, sim, pos_label=1)
            for ind2 in range(top_n):
            # for ind2 in arg_ind:
                if ri !=1:
                    i = sim>= sim[ind2]
                    j = sim <sim[ind2]
                #l2
                else:
                    i = sim <= sim[ind2]
                    j = sim > sim[ind2]
                #end l2
                tp = sum(labels[i])
                fn = sum(labels[j])
                tn = sum(labels[j]==0)
                fp = sum(labels[i]==0)
                # print(tp, fn)
                x.extend(fp/(fp+tn))
                y.extend(tp/k_dup) #tp / tp+fn
            if len(xs)==0:
                xs = np.array(x)
            else:
                xs +=x
            if len(ys) == 0:
                ys = np.array(y)
            else:
                ys += y

            # print(y)
        val = metrics.auc(x, y)
        print('val = ', round(val, 3))
        shift += 1 #下一个3重复组

    xs = xs/triples_num/k_dup
    ys = ys/triples_num/k_dup
    # print(xs[-1],ys[-1])
    # assert 1==0
    return xs, ys


def roc_curve(dataset, ri,pq, k_dup):
    xs = []
    ys = []
    #val_sample_ind 没有用
    # dataset = np.array(dataset, dtype='int64')

    dataset = np.array(dataset)

    dim = dataset.shape[1]
    data_num = dataset.shape[0]
    shift = 0
    triples_num = data_num//k_dup

    for i in range(triples_num):

        for ind in range(k_dup):
            # query_q = np.zeros((1,27))
            # query_q = np.zeros(dim)


            if ind != 0:
                continue
            query_q = dataset[ind+shift*k_dup].copy()
            # print(query_q)
            labels = np.zeros((data_num, 1))

            labels[shift * k_dup:shift * k_dup + k_dup] = 1
            # if ind == 0:
            #     labels[shift*k_dup+ind:shift*k_dup+ind+k_dup] =1
            # if ind ==1:
            #     labels[shift*k_dup+ind-1:shift*k_dup+ind+k_dup-1] = 1
            #     # labels[shift * 4 + ind] = 1
            #     # labels[shift * 4 + ind-1] = 1
            #     # labels[shift * 4 + ind+1] = 1
            #     # labels[shift * 4 + ind + 2] = 1
            # if ind ==2:
            #     labels[shift * k_dup + ind - 2:shift * k_dup + ind + k_dup - 2] = 1
            #     # labels[shift * 4 + ind + 1] = 1
            #     # labels[shift * 4 + ind] = 1
            #     # labels[shift * 4 + ind-1] = 1
            #     # labels[shift * 4 + ind-2] = 1
            # if ind == 3:
            #     labels[shift * k_dup + ind - 3:shift * k_dup + ind + k_dup - 3] = 1
            #     # labels[shift * 4 + ind - 3] = 1
            #     # labels[shift * 4 + ind] = 1
            #     # labels[shift * 4 + ind - 1] = 1
            #     # labels[shift * 4 + ind - 2] = 1
            # print('ind = ', ind + shift * 3, 'label = ', labels)
            if pq != None:
                print('pq:--')
                #asymmetric
                # dataset2 = pq.encode(dataset)
                # query_q = dataset[ind+shift*k_dup]

                #symmetric
                dataset2 = dataset.reshape(2*k_dup,-1)
                query_q = pq.decode(dataset[ind+shift*k_dup])[0]
                sim = linear_search2(dataset2, query_q, pq)
            else:
                res1,res2,sim ,topk_sim= linear_search(dataset,query_q,1,ri)
            sim = np.array(sim)
            arg_ind = np.argsort(sim)


            x=[0]
            y=[0]
            if ri != 1:
                arg_ind =arg_ind[::-1]

            # sim = sim[arg_ind[:6]]
            # labels = labels[arg_ind[:6]]
            # fpr, tpr, thresholds = metrics.roc_curve(labels, sim, pos_label=1)
            # for ind2 in range(len(arg_ind[:6])):
            for ind2 in arg_ind:
                if ri !=1:
                    i = sim>= sim[ind2]
                    j = sim <sim[ind2]
                #l2
                else:
                    i = sim <= sim[ind2]
                    j = sim > sim[ind2]
                #end l2
                tp = sum(labels[i])
                fn = sum(labels[j])
                tn = sum(labels[j]==0)
                fp = sum(labels[i]==0)
                # print(tp, fn)
                x.extend(fp/(fp+tn))
                y.extend(tp/(tp+fn)) #tp / tp+fn
            if len(xs)==0:
                xs = np.array(x)
            else:
                xs +=x
            if len(ys) == 0:
                ys = np.array(y)
            else:
                ys += y

            # print(y)
        shift += 1 #下一个3重复组

    xs = xs/triples_num/1
    ys = ys/triples_num/1
    # print(xs[-1],ys[-1])
    # assert 1==0
    return xs, ys


def concomitant_min_lsh(dataset, bits, dim, k,table_num):
    # n = 2**bits
    Ms = []
    bit = bits //table_num
    n1 = 2**bit

    tables = [{} for i in range(table_num)]

    for j in range(table_num):
        M = np.random.randn(dim, n1)
        Ms.append(M)
        for i in range(len(dataset)):
            P = np.dot(M.T, dataset[i])
            # min_ind = np.argsort(P)
            # min_inds = min_ind[:2]
            # for min_ind in min_inds:
            #     if min_ind not in tables[j]:
            #         tables[j][min_ind] = [i]
            #     else:
            #         tables[j][min_ind].append(i)
            # min_max
            min_ind = np.argsort(P)
            min_inds = min_ind[:k]
            max_inds = min_ind[-k:]
            # for min_ind in min_inds:
            #     min_ind = str(min_ind)+'min'
            #     if min_ind not in tables[j]:
            #         tables[j][min_ind] = [i]
            #     else:
            #         tables[j][min_ind].append(i)
            # for max_ind in max_inds:
            #     max_ind = str(max_ind) + 'max'
            #     if max_ind not in tables[j]:
            #         tables[j][max_ind] = [i]
            #     else:
            #         tables[j][max_ind].append(i)
            for mk in range(k):
                min_ind = min_inds[mk]
                max_ind = max_inds[mk]

                k4inds = []
                # k4ind = str(min_ind) + ',' + str(max_ind)
                # k4inds.append(k4ind)
                k4inds.append('min'+str(min_ind))
                k4inds.append('max'+str(max_ind))
                for k4i in k4inds:
                    if k4i not in tables[j]:
                        tables[j][k4i] = [i]
                    else:
                        tables[j][k4i].append(i)

    return tables,Ms

def len1Nor(v):
    return v/length(v)
def d2rank(data,flag=1):
    if flag ==1:
        for d in data:

            n0s = np.sum(d<=0)
            d2r = np.argsort(d)
            d[d2r] = range(1,38313+1)

            sum_rank_zeros = np.sum(range(1, n0s + 1))  # inmitate scipy.stats.stats.spearmanr

            d[d <= n0s] = sum_rank_zeros / n0s

            # d -= n0s  #original deal with zeros
            # d[d<0] = 0

        return data
    else:
        pass
def standardrized(d):
    mu = np.mean(d)
    sigmoid = np.std(d)
    # print(mu,sigmoid)
    return (d-mu)/sigmoid

def test(dataset, dataset0):    #some genes query
    q = [0 for i in range(38311)]
    q.extend([1,1])
    q.reverse()

    q = np.array(q)


    sims = []

    for d in dataset:
        sims.append(cos_sim(q, d))

    ind = np.argsort(sims)
    for i in ind[-10:]:
        print(sims[i],dataset[i][:2])
    assert 1==0

def normalization(dataset):
    data = []
    for d in dataset:
        data.append(d / np.linalg.norm(d))
    return np.array(data)



if __name__ == '__main__':

    import pandas as pd

    import scipy.signal
    mode= 'a'
    dataset,data_num, dim, sample_cond, sample_dup= read_data(mode)


    sam_cond = np.load('./sample_cond_10_21.npy')
    # dataset = np.load('./dataset_10_21.npy')
    # dataset = np.load('./data_rank_and_std_10_21.npy')
    # dataset = np.load('./data_lsh_104.npy',allow_pickle=True)
    # dataset = np.load('./dataset_pq_10_21.npy')
    # dels = []
    # i=0
    # for d in dataset:
    #     if np.sum(d>0.1) < 5000:
    #         dels.append(i)
    #         print(i)
    #     i+=1
    # dataset = np.delete(dataset, dels,axis=0)
    # sam_cond = np.delete(sam_cond,dels)






    sam_d = np.load('./sam_details_10_21.npy', allow_pickle=True).item()
    print(len(dataset), len(sam_cond))

    # i4sc = 0
    # for sc in sam_cond:
    #     if np.random.rand() < 0.1:
    #         print(dataset[i4sc][:10], sc)
    #     i4sc += 1
    # assert 1==0

    # sam_dup = {}
    # i = 0
    # for sc in sam_cond:
    #
    #     # print(sam_d[sc], sc)
    #     cd = sam_d[sc][1][0]
    #     gsm = sam_d[sc][0]
    #
    #     if cd[:-1] in sam_dup:
    #         ogsm = sam_dup[cd[:-1]][0][0]
    #         gsm = gsm.replace("GSM", '')
    #         ogsm = ogsm.replace('GSM', '')
    #         if abs(int(gsm)-int(ogsm)) > 100:
    #             continue
    #         sam_dup[cd[:-1]].append([sam_d[sc][0], cd, i])  #(gsm, title, ind)
    #     else:
    #         sam_dup[cd[:-1]] = [[sam_d[sc][0], cd, i]]
    #     i += 1
    #
    #
    #
    # np.save('./sam_dup_10_21',sam_dup)
    # assert 1==0

    #



    pq = None
    print(dataset.shape)

    # pq = nanopq.PQ(M=27,Ks=16)
    # pq.verbose = False
    # #
    # pq.fit(dataset)
    # # # # # # # #
    # with open('pq_27*4.pkl', 'wb') as f:
    #     pickle.dump(pq, f)
    # assert 1==0

    # # #
    # # #
    with open('pq_27*4.pkl', 'rb') as f:
        pq= pickle.load(f)  # pq_dumped is identical to pq
        pq.verbose = False


    #
    # #
    # dataset2 = []
    #
    # print(dataset.shape)
    # for d in dataset:
    #     dataset2.append(pq.encode(d.reshape(-1,38313)))
    # dataset2 = np.array(dataset2)
    # print(dataset2[0][:10])
    # np.save('./dataset_pq_10_21', dataset2)
    # assert 1==0
    # dataset2 = np.load('./data_lsh.npy', allow_pickle=True)
    # dataset3 = np.load('./dpq405.npy')

    sample_dup = np.load('./sam_dup_10_21.npy', allow_pickle=True).item()
    bits =104
    table_num = 1
    print(dataset.shape)

    #
    # print(tables)
    # assert 1==0
    # np.save('./tables', tables)
    # np.save('./planes', planes)
    # assert 100==1

    skip = False
    val_sample_ind = []
    k_dup = 7
    for k, v in sample_dup.items():    #v = [[gsm, title, ind],[gsm, ttl, ind],...]
        if len(v) >= k_dup:
            precd = ''
            tmplist = []
            repeat_gsm = {}
            for gsm,cd,ind in v:

                if len(precd) == 0:
                    precd = cd
                    tmplist.append(ind)
                    repeat_gsm[gsm] = 1
                    continue
                if cd[:-1] != cd[:-1] or  gsm in repeat_gsm:
                    print(v)    #not duplicate
                    break
                else:
                    tmplist.append(ind)
                    repeat_gsm[gsm] = 1
            else:

                val_sample_ind.append(tmplist)
    print(val_sample_ind)
    print(len(val_sample_ind))
    # assert 1==0

    # val_sample_ind = val_sample_ind[:50]
    # new_valsind = []
    # for vind in val_sample_ind:
    #     if np.sum(np.array(vind[:3])<2000) > 0:
    #         continue
    #     new_valsind.append(vind)
    # val_sample_ind = new_valsind



    xs = []
    ys = []

    cnt = len(val_sample_ind) * (len(val_sample_ind) - 1) / 2
    # --------------下面是roc曲线
    cnt = 0
    for i in range(len(val_sample_ind)):
        for j in range(i+1, len(val_sample_ind)):
            val_inds = []
            val_inds.extend(val_sample_ind[i][:k_dup])
            val_inds.extend(val_sample_ind[j][:k_dup])

            val_dataset = []
            dataset_pq = []
            dim = 38313
            for vi in val_inds[:]:
                # if pq !=None:
                #     val_dataset.append(pq.encode(dataset[vi].reshape(1,dim).astype('float32')))
                # else:
                val_dataset.append(dataset[vi])

            #
            #
            # val_dataset = d2rank(val_dataset)
            # for d in val_dataset:
            #     d[:] = standardrized(d)


            x, y = roc_curve(val_dataset, 1,pq, k_dup) #0--pearson, 1--l2, 2--spearman, 3--cos
            # x, y = roc_curve2(dataset[:500], val_inds, 1, pq)
            # print(type(x), x)
            val = metrics.auc(x, y)
            print(val)
            cnt+=1

            if len(xs) == 0:
                xs = x
            else:
                xs +=x
            if len(ys) == 0:
                ys = y
            else:
                ys +=y

    # #_________
    # val_inds = []
    # for i in range(len(val_sample_ind)):
    #     val_inds.extend(val_sample_ind[i][:k_dup])
    # val_dataset = []
    # for vi in val_inds:
    #     val_dataset.append(dataset[vi])
    # cnt =1
    # xs, ys = roc_curve2(val_dataset, 1, pq, k_dup)  # 0--pearson, 1--l2, 2--spearman, 3--cos

    # _________

    xs = xs / cnt
    ys = ys / cnt
    val = metrics.auc(xs, ys)
    np.save('x5', xs)
    np.save('y5', ys)
    print('avg val = ',val)
    strval = str(round(val, 3))
    plt.plot(xs, ys, '-o')
    plt.legend(['l2: ' + strval])
    plt.show()
    assert 1==0
    #--------------------------------
    # # ind4q = np.random.randint(0, data_num)


    # sample_cond = np.delete(sample_cond, ind4q)

    # print(tables)

    # for i in range(data_num):
    #     num_of_zeros = np.sum(dataset[i][:]==np.zeros((1,dim)))
    #     print(dim - num_of_zeros)

    min_k = 1
    Cascaded = 2  # cascaded = table
    top_k = 2
    dim = 38313
    tables, planes = lsh_tables(dataset, bits, dim, table_num)
    assert 1==0
    # tables2, Ms = concomitant_min_lsh(dataset, bits, dim, min_k, Cascaded)
    # np.save('./tables', tables)
    # np.save('./planes', planes)

    # tables = np.load('./tables.npy', allow_pickle=True)
    # planes = np.load('./planes.npy')

    # dataset = val_to_rank1(dataset)

    cnt = 0
    correct_cnt = 0
    total_time =0
    total_time2= 0
    qnum = len(dataset)//20
    # qinds = [np.random.randint(0, len(dataset)) for i in range(qnum)]
    qinds = [i for i in range(qnum)]
    # qinds = [3]
    for qind in qinds:
        # ind1 = np.random.randint(0,data_num)
        # ind2 = np.random.randint(0, data_num)
        # # print(ind1,ind2)
        #
        # query_q =dataset[ind1].copy()
        # # print('query_q:', ind1 ,sample_id[ind1], sample_cond[sample_id[ind1]])
        #
        # index = [i%2==0 for i in range(dim)]
        # query_q[index] = dataset[ind2][index]

        query_q = dataset[qind]
        # query_q0 = dataset0[qind]
        s1 = time.clock()
        total = []

        k = top_k
        topk_sim = [(0.0, -1) for i in range(k)]
        heapq.heapify(topk_sim)
        for i in range(table_num):

            sig4q = signature_bit(query_q, planes[i])

            bins = []
            bit4q = bin(sig4q)
            if len(bit4q)-2 < bits:
                zeors = ''.join(['0' for i in range(bits - len(bit4q)+2)])
                bit4q = zeors+bit4q[2:]
            else:
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
                    if cosine > topk_sim[0][0] and np.sum([p == ti[1] for ti in topk_sim]) == 0:
                        # print('before: ',topk_sim, cosine)
                        heapq.heappop(topk_sim)
                        heapq.heappush(topk_sim, (cosine, p))
                        # print('after:', topk_sim)
                    if cosine > max_cos:
                        max_cos = cosine
                        res = p
                # print(max_cos, res[:10])
                total.append((max_cos, res))

        if len(total) == 0:
            cnt +=1
            print('not found-----------------')
            continue
        res1 = max(total)
        # res1 = np.argsort(np.array(total)[:,0])
        e1 = time.clock()

        id = int(topk_sim[0][1])
        # print(id, correct_inds[cnt])
        s2 = time.clock()
        res2, cos2,total_sim,topk_l = linear_search(dataset, query_q, k)
        e2 = time.clock()

        if id != -1 and id == topk_l[0][1]:
            correct_cnt +=1

        condition = sample_cond[id]
        print(topk_sim[0], e1-s1,topk_l[0] , e2-s2)

        total_time += e1-s1
        total_time2 +=e2-s2
        cnt +=1
        # s2 = time.clock()
        # res2, cos2,total_sim = linear_search(dataset, query_q)
        # e2 = time.clock()
        #
        # id = int(res2)
        #
        # condition = sample_cond[id]
        # print(cos2, res2, e2 - s2, 'cond : '+condition)
        # print(dataset[res2][:10])
        # # print((e1-s1)/(e2-s2))
        # if res1[1] == res2:
        #     cnt += 1
        #
        # time1 +=e1-s1
        # time2 +=e2-s2

    print('presicion = ',correct_cnt/cnt, total_time/total_time2)
    # assert 1==0
    # 下面是Concomitant min hash algorithm
    #
    # cnt = 0
    # correct_cnt = 0
    # total_time = 0
    # total_time2=0
    #
    #
    # for qind in qinds:
    #     query_q = dataset[qind]
    #
    #
    #     s1 = time.clock()
    #     total_res = []
    #     k = top_k
    #     topk_sim = [(0.0, -1) for i in range(k)]
    #     heapq.heapify(topk_sim)
    #
    #     for i in range(Cascaded):
    #         # print(ind_min, tables2[i])
    #         M = Ms[i]
    #         p_of_q = np.dot(M.T, query_q)
    #         ind_min = np.argsort(p_of_q)
    #         ind_mins = ind_min[:min_k]
    #         ind_maxs = ind_min[-min_k:]
    #         # key4min = str(ind_min)+'min'
    #         # key4max = str(ind_max)+'max'
    #
    #
    #         for mk in range(min_k):
    #             k4inds = []
    #             # k4i = str(ind_mins[mk])+','+str(ind_maxs[mk])
    #             # k4inds.append(k4i)
    #             k4inds.append('min'+str(ind_mins[mk]))
    #             k4inds.append('max'+str(ind_maxs[mk]))
    #             for k4ind in k4inds:
    #                 if k4ind not in tables2[i]:
    #                     continue
    #                 max_cos = cos_sim(query_q, dataset[tables2[i][k4ind][0]])
    #                 res = 0
    #                 for ind in tables2[i][k4ind]:
    #                     cos = cos_sim(query_q, dataset[ind])
    #                     if cos > topk_sim[0][0] and np.sum([ind == i[1] for i in topk_sim]) == 0:
    #                         # print('before: ',topk_sim)
    #                         heapq.heappop(topk_sim)
    #                         heapq.heappush(topk_sim, (cos, ind))
    #                     if cos > max_cos:
    #                         max_cos = cos
    #                         res = ind
    #
    #         total_res.append((max_cos, res))
    #
    #
    #     if len(total_res) == 0:
    #         cnt += 1
    #         print('not found-----------------')
    #         continue
    #
    #
    #     res1 = max(total_res)
    #     s2 = time.clock()
    #
    #     id = int(topk_sim[0][1])
    #     # print(id, correct_inds[cnt])
    #     s3 = time.clock()
    #     res2, cos2, total_sim,topk_l = linear_search(dataset, query_q, k)
    #     s4 = time.clock()
    #     if id != -1 and id ==topk_l[0][1]:
    #         correct_cnt += 1
    #
    #
    #     # print(topk_sim[0],s2 - s1, topk_l[0], s4-s3)
    #     total_time +=s2-s1
    #     total_time2 +=s4-s3
    #     cnt+=1
    # print('presicion2 = ', correct_cnt / cnt, total_time/total_time2)