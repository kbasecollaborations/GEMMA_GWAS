import csv

with open('gemma.txt','r') as tsvin, open('gemma.tsv', 'w') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    tsvout = csvout.write('SNP\tCHR\tBP\tP\n')

    c = 0
    for row in tsvin:
        if c != 0 and int(row[0]) == 1:
            # row[0] chr
            # row[1] SNP
            # row[2] pos
            # row[8] pval:
            csvout.write(row[1]+'\t'+row[0]+'\t'+row[2]+'\t'+row[8]+'\n')
        c = c+1
