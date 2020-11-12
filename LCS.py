import numpy as np
def seq_align(x, y):
    m = len(x); n = len(y)
    if (m == 0):
        return d * n
    if (n == 0):
        return d * m

    if (m,n) in dic:
        return dic[(m, n)]

    score1 = seq_align(x[:m-1], y[:n-1]) + int(x[-1]==y[-1])
    score2 = seq_align(x[:m-1], y) + d
    score3 = seq_align(x, y[:n-1]) + d

    max_score = max(max(score1,score2), score3)
    if max_score == score1:
        b[m-1][n-1] = 7
    elif max_score == score2:
        b[m-1][n-1] = 8
    else:
        b[m-1][n-1] = 4

    dic[(m,n)]= max_score
    return max_score

def seq_print(b, x, y, s1, s2):
    m = len(x)
    n = len(y)
    if (m == 0):
        if (n > 0):
            s1 += '-'*n
            s2 += y[0:n]
        return

    if (n == 0):
        if (m > 0):
            s1 += x[0:m]
            s2 += '-'*m
        return

    if b[m-1][n-1] == 7:
        seq_print(b, x[:m-1], y[:n-1], s1, s2)
        s1 += x[-1]
        s2 += y[-1]

    elif b[m-1][n-1] == 4:
        seq_print(b, x, y[:n-1], s1, s2)
        s1 += '-'
        s2 += y[-1]

    else:
        seq_print(b, x[:m-1], y, s1, s2)
        s1 += x[-1]
        s2 += '-'

if __name__ == '__main__':
    # a = ['A', 'C', 'G', 'T']

    # L = [[2, -7, -5, -7], \
    #      [-7, 2, -7, -5], \
    #      [-5, -7, 2, -7], \
    #      [-7, -5, -7, 2]]

    # linear gap penalty
    d = 0

    #dic存储长度分别为m,n的字符串的得分，(m,n)为key， 得分为value


    xs =['ColandINA-RPKM', 'Ht1-5andINA-RPKM', 'ColandINA-RPKM', 'Nr1andINA-RPKM', 'Col-RPKM', 'ColandINA-RPKM']
    ys = ['Col', 'hac1/5', 'npr1-1', 'hac1/5+INA', 'Col+INA', 'npr1-1+INA']
    for x in xs:
        print('---')
        for y in ys:
            #b[i,j]记录长度(i,j)的得分的来源。b用于重构最优解
            b = [[5 for i in range(len(y))] for i in range(len(x))]
            dic = {}
            seq_align(x, y)
            s1 = []
            s2 = []
            seq_print(b, x, y, s1, s2)

            str1 =''
            for s in s1:
                str1 += s
            str2 =''
            for s in s2:
                str2 += s

            # print(str1)
            # print(str2)
            print(np.sum([int(str1[i]==str2[i] and str1[i]!='-') for i in range(len(str1))]), end='\t')
            #print(b)

