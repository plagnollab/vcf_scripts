#! /bin/env python
from __future__ import print_function
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-m','--mapping.csv', dest="mapping", help="CSV mapping file", metavar="")
(options, args) = parser.parse_args()

m=options.mapping 
mapping=dict([l.strip().split(',') for l in file(m,'r').readlines()])

for l in sys.stdin:
    l=l.strip()
    # skip lines which start with ##
    #ignore all other '##' lines
    if l.startswith('##'):
        print(l)
        continue
    #line which starts with a single '#' is the header
    #When we get to the header line, print it out.
    if l.startswith('#'):
        l=l.strip('#')
        HEADER=l.strip().split('\t')
        HEADER=[ mapping[h] if h in mapping else h for h in HEADER ]
        print('#','\t'.join(HEADER))
        continue
    print(l)



