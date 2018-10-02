#Jensen-Shannon Divergence 
#Created by Qian Zhang 07/16/2014
#revised on 12/30/2014 

#!/usr/bin/python


import sys, warnings, os

sys.path.append('/shared/apps/python/python-2.7.6/lib/python2.7/')
sys.path.append('/shared/apps/python/python-2.7.9/lib/python2.7/site-packages')
import numpy as np 



#JS Divergence calculation function jsd()#
def jsd(fileone,filetwo):
    kmerone = []
    valueone = []
    kmertwo = []
    valuetwo = []
    kmer = []
    fp = []
    fq = []
    p = []
    q = []
    with open (fileone, 'r') as infile1:
        for line in infile1:
            kmerone.append(line.strip().split(' ')[0])
            valueone.append(line.strip().split(' ')[1]) 
    with open (filetwo, 'r') as infile2:
        for item in infile2:
            kmertwo.append(item.strip().split(' ')[0])
            valuetwo.append(item.strip().split(' ')[1])
    kmer = list(set(kmerone + kmertwo))
    dict1 = dict(zip(kmerone,valueone))
    dict2 = dict(zip(kmertwo,valuetwo))
    for k in kmer:
        if dict1.has_key(k):
            fp.append(int(dict1[k]))
        else:
            fp.append(0)
    for k in kmer:
        if dict2.has_key(k):
            fq.append(int(dict2[k]))
        else:
            fq.append(0)
    p = list(np.array(fp)/(float(sum(fp))))
    q = list(np.array(fq)/(float(sum(fq))))
    m=(np.array(p)+np.array(q))/2
    value = (0.5 * sum(_p * np.log(_p / _m) for _p, _m in zip(p, m) if _p != 0) + \
             0.5 * sum(_q * np.log(_q / _m) for _q, _m in zip(q, m) if _q != 0))/np.log(2)
    return value



## =================================================================
## argument parser

import argparse

parser = argparse.ArgumentParser(description="JS Divergence with python",
    prog = 'python', #program name
    prefix_chars='-', # prefix for options
    fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
    conflict_handler='resolve', # for handling conflict options
    add_help=False, # include help in the options
    formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
    )

## input bc files and/or directories that include bc files to be pairwise calculate#

parser.add_argument("-f",help="first file",dest='firstfile')
parser.add_argument("-s",help="second file",dest='secondfile')
parser.add_argument("-o",help="output file",dest='outputfile')
## ===========
## Main function
## ===========
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)
    with open (args.outputfile, 'w')as output:
        output.write(args.firstfile.strip('.fna.txt') +'\t' + args.secondfile.strip('.fna.txt') +'\t' + str(jsd(args.firstfile, args.secondfile)))

    

#run
if __name__ == "__main__":
    sys.exit(main())

