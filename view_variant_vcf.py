import sys

#View single variant identified by rsID:
#python view_variant.py file.txt r [rsID]

#View single variant identified by chromosome and position:
#python view_variant.py file.txt p [chr] [pos]

input=open(sys.argv[1],"r")
lookup=str(sys.argv[2])

if lookup=="r":
	rs=str(sys.argv[3])
if lookup=="p":
	chromosome=sys.argv[3]
	pos=sys.argv[4]

array=[]

first=0
for line in input:
	split=line.rstrip("\n").split("\t")
	if first==1:
		if lookup=="r":
			if split[2]==rs:
				array.append(split)
		if lookup=="p":
			if split[0]==chromosome and split[1]==pos:
				array.append(split)
        if split[0]=="#CHROM":
                numbers=[]
                for i in range(len(split)):
                        numbers.append(str(i+1))
                array.append(numbers)
                array.append(split)
		first=1

transposed=zip(*array)

for line in transposed:
	print "\t".join(line)
