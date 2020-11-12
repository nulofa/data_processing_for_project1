exp1 =[]
exp2 = []
from scipy.stats.stats import spearmanr,pearsonr
with open('./sampel_exp.results') as f:
    for line in f:
        cont = line.split('\t')
        if cont[0].startswith('gene-'):
            gn = cont[0].split('-')[1].strip()
            exp1.append(float(cont[-2]))
            exp2.append(float(cont[-1]))


print(spearmanr(exp1, exp2))
print(pearsonr(exp2, exp1))
