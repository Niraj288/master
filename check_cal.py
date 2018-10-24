#!/usr/bin/python
import os
import module
import sys

def check_error(path):
	f=open(path[:-4]+'.error')
	lines=f.readlines()
	f.close()
	if 'CANCELLED' in ''.join(lines):
		return 'TIME LIMIT EXCEEDED!! in '
	return 'Probably running/terminated for '

def func(path):
        #print path 
        ref,ref_f=0,1
        file=open(path,'r')
        file_=file.readlines()
        file.close()
    
	if 'Normal termination of Gaussian' in file_[-1] and '-done' in sys.argv:
                print 'Normal termination in '+path
	if '-c' in sys.argv and 'Convergence failure' in ''.join(file_[-50:]):
                print 'Convergence failure in '+path
	if 'ORCA TERMINATED NORMALLY' in file_[-2] and '-done' in sys.argv:
                print 'Normal termination in '+path
	if 'O   R   C   A ' in ''.join(file_[:1000]) and 'ORCA TERMINATED NORMALLY' not in file_[-2]:
		print 'Error termination in '+path
		return
	if 'Gaussian' not in ''.join(file_[0:5]):
                return
	if '-done' in sys.argv or '-c' in sys.argv:
		return
        for line in file_:
                if 'Frequencies' in line:
                        k=line.strip().split()
                        if len(k)>3:
				if float(k[3])<0:
                                	ref_f=2
                                	break
                        if float(k[2])<0:
                                ref_f=0
                                break
	c_ref=1
        if ref_f==0:
                print '1 Imaginiary frequencies in '+path
        elif ref_f==2:
                print 'More than 1 Imaginiary frequencies in '+path
        elif 'Convergence failure' in ''.join(file_[-50:]):
		print 'Convergence failure in '+path
		c_ref=0
	if 'Normal termination of Gaussian' not in file_[-1]:
		if c_ref:
			print check_error(path)+path
                ref_f=3
        if ref_f!=1:
                return 0
        return 1

if __name__=='__main__':
        p=''
        if len(sys.argv)>1:
        	p=sys.argv[1]
        path=p or raw_input("Enter path : ")
        module.search_deep(path,func,['out'])




