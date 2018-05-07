#!/usr/bin/python

import os
import sys
import threading
import multiprocessing
d={'ss':'standard-mem-s',
'sm':'standard-mem-m',
'sl':'standard-mem-l',
'ms':'medium-mem-1-s',
'mm':'medium-mem-1-m',
'ml':'medium-mem-1-l',
'm2':'medium-mem-2',
'h1':'high-mem-1',
'h2':'high-mem-2',
'gpu1':'gpgpu-1',
'gpu2':'gpgpu-2',
'mic':'mic',
'htc':'htc'}

def gaussian(path,d,filen):
	filename,en='.'.join(filen.strip().split('.')[:-1]),filen.strip().split('.')[-1]
        os.system('mkdir '+filename)
        os.chdir(filename)
        os.system("cp ../"+filen+' .')
	jojo=open(filename+'.job','w')
	s1="""#!/bin/bash

#SBATCH -J """
	jn=filename
	s12="""
#SBATCH -o nn_out.out
#SBATCH -p """
	s2=d[sys.argv[2]]
	s3="""
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --mem=100000
"""
	if en=='rub':
		s4="python /users/nirajv/side_projects/rub_data_ml.py "+filen
	elif en=='kat':
		s4="python /users/nirajv/side_projects/katja_ml.py "+filen
	else:
		print 'Input not recognized!!'
		return
	jojo.write(s1+jn+s12+s2+s3+s4)
	jojo.close()
        os.system('sbatch '+filename+'.job')


gaussian('',d,sys.argv[1])
