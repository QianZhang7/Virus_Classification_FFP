import sys
sys.path.append('/shared/apps/python/python-2.7.6/lib/python2.7/')
import numpy as np

#Relative Entropy Calculation Function rec()#
def rec(folder, filename, lower, upper):
    files = []
    for i in range((int(lower)-2), (int(upper) +1)):
        files.append(folder + "kmer"+ str(i) +"/" + filename)
    #initialize values
    value = []
    countsum = []
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
        if x >=2:
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
            re = sum(_p* np.log(_p / _m) for _p, _m in zip(fi, fihat))/np.log(2)
            value.append(re)
    return value


## =================================================================
## argument parser

import argparse

parser = argparse.ArgumentParser(description="Relative entropy with python",
        prog = 'python', #program name
        prefix_chars='-', # prefix for options
        fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
        conflict_handler='resolve', # for handling conflict options
        add_help=False, # include help in the options
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
        )

## input bc files and/or directories that include bc files to be pairwise calculate#

parser.add_argument("-d",help="directory",dest='directory')
parser.add_argument("-f",help="filename",dest='filename')
parser.add_argument("-o",help="output file",dest='outputfile')
parser.add_argument("-l",help="lowerK",dest='lowerK')
parser.add_argument("-u",help="upperK",dest='upperK')

## ===========
## Main function
## ===========
def main(argv=None):
    result = []
    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)
    with open (args.outputfile, 'w')as output:
        output.write (args.filename.strip('.fna.txt') + '\t')
        result.append(rec(args.directory, args.filename, args.lowerK, args.upperK))
        for r in result:
            output.write("%s\n" % r)

#run
if __name__ == "__main__":
    sys.exit(main())
