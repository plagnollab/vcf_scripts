#! /bin/env python
from __future__ import print_function
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-s','--sample', dest="sample", action='append', default=[], help="", metavar="")
parser.add_option('-v','--variant', dest="variant", action='append', default=[], help="", metavar="")
(options, args) = parser.parse_args()


for l in sys.stdin:
    l=l.strip()
    # skip lines which start with ##
    #ignore all other '##' lines
    if l.startswith('##'): continue
    #line which starts with a single '#' is the header
    if l.startswith('#'):
        l=l.strip('#')
        HEADER=l.strip().split('\t')
        SAMPLES=HEADER[9:]
        GENOTYPE_HEADER=['VARIANT_ID']+SAMPLES
        continue
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))
    for k in SAMPLES: s[k]=s[k]
    VARIANT_ID='_'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']])
    s['VARIANT_ID']=VARIANT_ID
    if VARIANT_ID in options.variant:
        for samp in options.sample:
            print(VARIANT_ID, samp, s[samp], sep=' ')




