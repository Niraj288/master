import os
import module

def func(path):
	ref,ref_f=0,1
	file=open(path,'r')
	file_=file.readlines()
	file.close()
	for line in file_:
		if 'Frequencies' in line:
			k=line.strip().split()
			if float(k[2])<0:
				ref_f=0
				break
	if ref_f==0:
		print 'Imaginiary frequencies in '+path
	if 'Normal termination of Gaussian' not in file_[-1]:
		print 'Error termination in '+path
	

path=raw_input("Enter path : ")
module.search_deep(path,func,['g16.out','g09.out'])