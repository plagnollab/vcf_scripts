### Hopefully vcfbreakmulti, from https://github.com/ekg/vcflib, can take care of this now.
#! /bin/env python
from __future__ import print_function
import sys
import argparse 
import os.path
from mmap import mmap
import time
from multiprocessing import Process
import subprocess
import csv

usage=""" """

#these 9 column headers are standard to all VCF files
STD_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
#the samples headers depend on the number of samples in the file
#which we find out once we read the #CHROM line
SAMPLE_HEADERS=[]

parser=argparse.ArgumentParser(description='Arguments to vcftools.py: this script will be default skip variants which are missing in all samples.')
parser.add_argument('--GQ',default=None,type=int,help='')
parser.add_argument('--DP',default=None,type=int,help='set variants with depth < threshold to NA')
parser.add_argument('--samples',default=None,type=str,help='file containing samples')
parser.add_argument('--vcf',default=False,action='store_true',help='VCF output with single allele variants, subsets if sample specified')
parser.add_argument('--freq',default=False,action='store_true',help='freq output, subsets if sample is specified')
parser.add_argument('--pileup',default=False,action='store_true',help='pileup output, only output variants')
args=parser.parse_args()

if not args.vcf and not args.freq and not args.pileup: 
    print(help())
    exit(1)

if args.samples: SAMPLE_HEADERS=[l for l in csv.DictReader(file(args.samples,'r'))]

def print_line(s):
    print( *([s[h] for h in STD_HEADERS+SAMPLE_HEADERS]), sep='\t' )

def frequency(s):
    geno=[s[h].split(':')[0] for h in SAMPLE_HEADERS]
    total_count=len(geno)
    missing_count=geno.count('./.')
    wt_count=geno.count('0/0')
    het_count=geno.count('0/1')+geno.count('1/0')
    hom_count=geno.count('1/1')
    N=total_count-missing_count
    if N==0: AF=0
    else: AF=float(2*hom_count+het_count)/(2*N)
    return {'N':N, 'total_count':total_count, 'missing_count':missing_count, 'wt_count':wt_count, 'het_count':het_count, 'hom_count':hom_count, 'AF':AF }

def print_frequency(s):
    f=frequency(s)
    N=f['N']
    AF=f['AF']
    hom_count=f['hom_count']
    het_count=f['het_count']
    wt_count=f['wt_count']
    total_count=f['total_count']
    missing_count=f['missing_count']
    VARIANT_ID='_'.join([s['CHROM'], s['POS'], s['REF'], s['ALT']])
    if N!=0: print(VARIANT_ID,total_count,missing_count,wt_count,het_count,hom_count,AF, sep=',' )

def print_pileup(s):
    f=frequency(s)
    N=f['N']
    if N!=0: print(s['CHROM'], s['POS'], s['REF'], s['ALT'], sep=' ')

if args.freq:
    print('VARIANT_ID','total_count','missing_count','wt_count','het_count','hom_count','AF',sep=',')
elif args.vcf: 
    print('##fileformat=VCFv4.1')

for line in sys.stdin:
    line=line.strip()
    #remove beginning
    #lines which start with '###'
    #are not tab separated
    if line.startswith('##'):
        if args.vcf: print(line)
        continue
    #this is tab separated line
    s=line.split("\t")
    #header line yay!: #CHROM ...
    if line.startswith('#'):
        if args.vcf: print(line)
        headers=s
        headers[0]=headers[0].strip('#')
        #the first 9 names in the header are standard (see above)
        if (headers[0:len(STD_HEADERS)] != STD_HEADERS): raise 'hell'
        #everything else in the header is a sample name
        if not len(SAMPLE_HEADERS): SAMPLE_HEADERS=headers[len(STD_HEADERS):]
        continue
    #you can now access elements by their header name
    s=dict(zip(headers,s))
    #I would expect each sample to be formatted according to the format column.
    #However I found in practise this not true.
    #Only the first two fields: GT (genotype) and AD (allele depth)
    #are compulsory, the remaining fields, including DP (total depth),
    #can be missing.
    #So I calculate DP from AD instead.
    #I ignore the other fields for now.
    for h in SAMPLE_HEADERS:
        d=dict(zip(s['FORMAT'].split(':'),s[h].split(':')))
        #print(s['FORMAT'])
        #print(d)
        GT=d['GT']
        AD=d['AD']
        #DP=d['DP']
        DP=str(sum([int(x) for x in AD.split(',')]))
        #GQ can be missing
        GQ=d.get('GQ','.')
        # if GQ is '.' the must be missing
        if GQ == '.' and GT != './.': raise 'hell'
        # if fails QC set to missing
        if ((args.DP and int(DP) < args.DP) or (args.GQ and GQ!='.' and int(GQ) < args.GQ)): GT='./.'
        PL=d.get('PL','.')
        s[h]=':'.join([GT,AD,DP,GQ,PL])
    #split alternate alleles
    alternative_alleles=s['ALT'].split(',')
    #if single alternate allele
    #then just print out as normal
    if len(alternative_alleles)<=1:
        s['FORMAT']='GT:AD:DP:GQ:PL'
        if args.freq: print_frequency(s)
        elif args.vcf: print_line(s)
        elif args.pileup: print_pileup(s)
        continue
    #otherwise split over as many lines as there are alternative alleles
    #each genotype then takes a 1 if matches the alternative or a 0 otherwise
    alleles=[s['REF']]+alternative_alleles
    n2geno=dict(zip(['.']+[str(i) for i in range(0,len(alleles))],['.']+alleles))
    for idx, alt in enumerate(alternative_alleles):
        s['ALT']=alt
        #recode GT
        s1=s.copy()
        for h in SAMPLE_HEADERS:
            d=dict(zip(s1['FORMAT'].split(':'),s1[h].split(":")))
            #length of allele depth is 2 where first is always REF allele depth
            #and second can be either ALT
            AD=d['AD'].split(',')
            AD[1]=AD[idx+1]
            AD=','.join(AD[:2])
            # matches REF -> 0
            # matches ALT -> 1
            # matches anything else -> .
            #sorted so that "." precedes number
            GT='/'.join(sorted([{s['REF']:'0',s['ALT']:'1'}.get(n2geno[g],'.') for g in d['GT'].split('/')]))
            DP=str(sum([int(x) for x in AD.split(',')]))
            s1[h]=':'.join( [GT, AD, DP] )
            s1['FORMAT']='GT:AD:DP'
        if args.freq: print_frequency(s1)
        elif args.vcf: print_line(s1)
        elif args.pileup: print_pileup(s1)





