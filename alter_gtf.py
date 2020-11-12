new = []
with open('./Species.genome.gtf') as f:
    for line in f:
        if 'gene_id' not in line.split('\t')[8]:
            newline = line.replace('\n', 'gene_id "null";\n')
            new.append(newline)
        else:
            new.append(line)

with open('./nnj.gtf', 'w') as f:
    for line in new:
        f.write(line)