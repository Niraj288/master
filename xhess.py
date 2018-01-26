#!/usr/bin/python
import os
import threading
import multiprocessing

def gaussian(path):
        filename=path
        #os.system('mkdir '+filename)
        os.chdir(filename)
        #os.system("cp ../"+path+' ./')
	os.system('sbatch xrunh -q gpgpu-1 -f xtbopt.xyz')
processes=[]
for i in os.listdir('.'):
        if os.path.isdir(i):
                p=multiprocessing.Process(target=gaussian, args=(i,))
                processes.append(p)

for p in processes:
        p.start()


print 'All calculations are done!!!'
