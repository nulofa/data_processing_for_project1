from __future__ import division
import numpy as np
import math
import time

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


def bitcount(n):
    """
    gets the number of bits set to 1
    """
    count = 0
    while n:
        count += 1
        n = n & (n - 1)
    return count

def length(v):
    """returns the length of a vector"""
    return math.sqrt(np.dot(v, v))

def cos_sim(v1, v2):
    assert len(v1) == len(v2)
    return np.dot(v1,v2)/ length(v1)/length(v2)

def linear_search(data, q):
    max_cos= cos_sim(data[0], q)
    res = 0
    data_num = len(data)
    for p in range(data_num):
        cosine = cos_sim(data[p], q)
        if cosine > max_cos:
            max_cos = cosine
            res = p
    return res, max_cos

def lsh_tables(data, bits, dim,iter):
    tables = []
    planes = []
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

if __name__ == '__main__':
    dim = 10000  # dimension of data points (# of features)
    bits = 9    # number of bits (planes) per signature

    data_num = 1000

    data = np.random.randn(data_num, dim)

    # print(data2[0])
    table_num = 20

    tables, planes = lsh_tables(data, bits, dim, table_num)

    string = input('go 4 start, end 4 end')
    total_res = []
    while string != 'end':
        query_num = 100
        cnt = 0
        dif = 0
        for qn in range(query_num):
            query_q = np.random.randn(dim)

            s1 = time.clock()
            total = []
            for i in range(table_num):

                sig4q = signature_bit(query_q, planes[i])

                bins = []
                bit4q = bin(sig4q)
                bit4q = bit4q[2:]
                bins.append(int(bit4q,2))
                # print('bit4q= ',bit4q)
                for b in range(len(bit4q)):
                    if bit4q[b] == '0':
                        # print(bit4q[:b]+'1'+bit4q[b+1:])
                        bins.append(int(bit4q[:b]+'1'+bit4q[b+1:], 2))
                    else:
                        bins.append(int(bit4q[:b] + '0' + bit4q[b + 1:], 2))


                # if sig4q not in tables[i]:
                #     continue
                # print('len = ', len(tables[i][sig4q]))

                for sig4q in bins:
                    if sig4q not in tables[i]:
                        continue
                    ind = tables[i][sig4q][0]
                    max_cos = cos_sim(query_q, data[ind])
                    res = ind

                    for p in tables[i][sig4q]:
                        cosine = cos_sim(data[p], query_q)
                        if cosine > max_cos:
                            max_cos = cosine
                            res = p
                    #print(max_cos, res[:10])
                    total.append((max_cos, res))


            res1 = max(total)
            e1 = time.clock()
            print(res1, e1-s1)

            s2 = time.clock()
            res2, cos2 = linear_search(data, query_q)
            e2 = time.clock()
            print(cos2, res2, e2-s2)
            dif += cos2 - res1[0]
            if res1[1] == res2:
                cnt +=1

        print(cnt/query_num, dif)
        total_res.append(cnt/query_num)
        string = input('go 4 start, end 4 end')

    print(np.average(total_res))
