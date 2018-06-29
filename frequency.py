import os
import sys
if len(sys.argv)>1:
	path=sys.argv[1]
else:
	path=raw_input("Enter path : ")
file =open(path,'r')

lis=file.readlines()
file.close()
l1,l2=[],[]
for i in range (len(lis)):
	if 'Frequencies' in lis[i]:
		l1+=lis[i].strip().split()[2:]
		l2+=lis[i-1].strip().split()


for i in range (len(l1)):
	print i+1,l1[i].split('.')[0],l2[i]
