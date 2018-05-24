#!/usr/bin/python
import os
import threading
import multiprocessing
import time
import sys
directories={}

def threa(command):
        os.system(command)

def gauss_claculations(path):
        refe=1
	if len(sys.argv)>1:
		keyword=sys.argv[1]
	else:
        	keyword='.wfx'
	processes=[]
        for i in os.listdir(path):
                if  keyword==i[-len(keyword):]:
                        if refe==0:
                                if raw_input('Submit job for '+i+' ? (y/n) : ')=='y':
                                        filename=i.strip().split('.')[0]
                                        #mode='mm'
                                        directories[filename]=1
                                                 #threading.Thread(target=threa, args=(path,mode,filename,)).start()
                                        p=multiprocessing.Process(target=threa, args=("/hpc/applications/aimall/17.01.25/aimqb.ish -nogui "+i,))
                                        processes.append(p)
                        else:
                                filename=i.strip().split('.')[0]
                                #mode='mm'
                                directories[filename]=1
                                threading.Thread(target=threa, args=("/hpc/applications/aimall/17.01.25/aimqb.ish -nogui "+i,)).start()
                                #p=multiprocessing.Process(target=threa, args=("aimqb.ish -nogui "+i,))
                                #processes.append(p)

        for p in processes:
                p.start()

        for p in processes:
                p.join()

        print 'All jobs Done !! \n'
path='.'
gauss_claculations(path)
