import sys, warnings, os
sys.path.append('/shared/apps/python/python-2.7.6/lib/python2.7/site-packages')
import pandas

kmerdict = {}
lower = 5
upper = 15
folder =  "/compgenpanfs/qzo/virusbig_jan2015/result_jan/kmer_data/"
filename = "AC_000004.fna.txt"


countsum = []
files = []
for i in range((lower-2), (upper +1)):
    files.append(folder + "kmer"+ str(i) +"/" + filename)

value = []
kmerdict = {}
for x, item in enumerate(files):
    kmer = []
    count = []
    with open(item, 'r') as filein:
        for line in filein:
            kmer.append(line.strip().split(' ')[0])
            count.append(int(line.strip().split(' ')[1]))
    countsum.append(sum(count))
    kmerdict[x]= dict(zip(kmer, (count)))
    if x >= 2:
        fi = []
        fihat = []
        ful = []
        fur = []
        fd = []
        re = []
        for k in kmerdict[x]:
            fi.append(int(kmerdict[x][k])/(float(countsum[x])))
            ful.append(int(kmerdict[(x-1)][k[1:]])/(float(countsum[x-1])))
            fur.append(int(kmerdict[(x-1)][k[:-1]])/(float(countsum[x-1])))
            fd.append(int(kmerdict[(x-2)][k[1:-1]])/(float(countsum[x-2])))
        fihat = list(numpy.array(ful)*numpy.array(fur)/(numpy.array(fd)))
        re = sum(_p* math.log(_p / _m) for _p, _m in zip(fi, fihat))/math.log(2)
        value.append(re)

print x, len(kmerdict), countsum, value

