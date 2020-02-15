#!/usr/bin/env python3

'''

The Cut & tag has a lot more Transposase sequences in, need to trim them out:

>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

'''

import sys, os
import gzip as gzipfile

def fastqPE(filename1, filename2, gzip=True):
    """
    generator object to parse fastQ PE files

    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((

    """
    if gzip:
        oh1 = gzipfile.open(filename1, "rt")
        oh2 = gzipfile.open(filename2, "rt")
    else:
        oh1 = open(filename1, "rt")
        oh2 = open(filename2, "rt")

    name1 = "dummy"
    while name1 != "":
        name1 = oh1.readline().strip()
        seq1 = oh1.readline().strip()
        strand1 = oh1.readline().strip()
        qual1 = oh1.readline().strip()

        name2 = oh2.readline().strip()
        seq2 = oh2.readline().strip()
        strand2 = oh2.readline().strip()
        qual2 = oh2.readline().strip()

        res = ({"name": name1, "strand": strand1, "seq": seq1, "qual": qual1},
            {"name": name2, "strand": strand2, "seq": seq2, "qual": qual2})
        yield res
    return

adapters = {
    'Trans1':     'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
    'Trans1_rc':  'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
    'Trans2':     'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
    'Trans2_rc':  'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
    }
counts = {k: 0 for k in adapters}
counts['rejected_reads'] = 0

adapters = {v: k for k,v in adapters.items()}

print(adapters)

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('trimit.py infileP1 infileP2 outfileP1 outfileP2')

    infile_p1 = sys.argv[1]
    infile_p2 = sys.argv[2]
    outfile_p1 = sys.argv[3]
    outfile_p2 = sys.argv[4]

    p1 = gzipfile.open(outfile_p1, "wt")
    p2 = gzipfile.open(outfile_p2, "wt")

    for index, read in enumerate(fastqPE(infile_p1, infile_p2)):
        seq1 = read[0]['seq']
        seq2 = read[1]['seq']
        qual1 = read[0]['qual']
        qual2 = read[1]['qual']
        for a in adapters:
            k = adapters[a]
            if a in seq1:
                p = seq1.find(a)
                seq1 = seq1[0:p]
                qual1 = qual1[0:p]
                counts[k] += 1

            if a in seq2:
                p = seq2.find(a)
                seq2 = seq2[0:p]
                qual2 = qual2[0:p]
                counts[k] += 1

        if len(seq1) < 50 and len(seq2) < 50:
            counts['rejected_reads'] += 1
        else:
            p1.write('{0}\n{1}\n+\n{2}\n'.format(read[0]['name'], seq1, qual1))
            p2.write('{0}\n{1}\n+\n{2}\n'.format(read[1]['name'], seq2, qual2))

        if (index+1) % 1e5 == 0:
            print('{:,}'.format(index+1))

    p1.close()
    p2.close()

for k in counts:
    print('{0}: {1:,}'.format(k, counts[k]))

