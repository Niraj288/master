#!/usr/bin/python
import os
import threading
import multiprocessing

def gaussian(path):
        filename=path[:-4]
        os.system('mkdir '+filename)
        os.chdir(filename)
        os.system("cp ../"+path+' ./')
        os.system('sbatch xrun -q gpgpu-1 -f '+path)
processes=[]
for i in os.listdir('.'):
        if i[-4:]=='.xyz':
                p=multiprocessing.Process(target=gaussian, args=(i,))
                processes.append(p)

for p in processes:
        p.start()


print 'All calculation are done!!!'
