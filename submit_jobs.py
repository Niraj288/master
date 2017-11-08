#!/usr/bin/python
#Automation
import os
import threading
import multiprocessing
import datetime
import sys
import time
directories={}
referance=0
mod=''
def schedule(st):
	global mod
	mod=sys.argv[sys.argv.index('-m')+1]
        h,m,s=map(int,st.split(':'))
        now = datetime.datetime.now()
        today = now.replace(hour=h, minute=m, second=s, microsecond=0)
        while now<today:
                time.sleep(3)
		now = datetime.datetime.now()
	su()
def su():
	global referance
	print 'submitting job at '+str(datetime.datetime.now())
	referance=1 
        gauss_claculations(os.listdir(os.getcwd()))

def threa(command):
        os.system(command)

def gauss_claculations(lis_paths):
	global directories,referance,mod
	path,filename='',''
	if referance==1:
		mode=mod
	else:
		if '-m' in sys.argv:
			ind=sys.argv.index('-m')
			mode=sys.argv[ind+1]
		else:
			mode=raw_input("Enter mode : ")
	processes=[]
	refe=0
	if '-all' in sys.argv or referance==1:
		refe=1
	elif raw_input('Send all .g16 or .inp files in the directory ?? (y/n) : ')=='y':
		refe=1
	if len(lis_paths)>0:
		refe=1
		paths=lis_paths
	else:
		paths=os.listdir(os.getcwd()) 
	for i in paths:
	        if '.g16'==i[-4:] or '.inp'==i[-4:]:
	        	if refe==0:
	        		if raw_input('Submit job for '+i+' ? (y/n) : ')=='y':
		               		 filename=i.strip().split('.')[0]
		               		 #mode='mm'
		               		 directories[filename]=1
		               		 path=filename+'.g16'
		               		 #threading.Thread(target=threa, args=(path,mode,filename,)).start()
		               		 p=multiprocessing.Process(target=threa, args=("grun "+path+' '+mode,))
		               		 processes.append(p)
		        else:
		        	filename=i.strip().split('.')[0]
		                #mode='mm'
		                directories[filename]=1
		                path=filename+'.g16'
		                #threading.Thread(target=threa, args=(path,mode,filename,)).start()
		                p=multiprocessing.Process(target=threa, args=("grun "+path+' '+mode,))
		                processes.append(p)

	for p in processes:
	        p.start()

	for p in processes:
	        p.join()

	print 'All jobs sent !! \n'

def organize():
	global directories
	l_files=os.listdir(os.getcwd())
	for file in directories:
		if file in l_files:
                        print file,"directory already present!!\n All files with the prefix "+file+" will be ignored!! "
		else:
			os.system("mkdir "+file)
	for file in l_files:
		f=file.split('.')[0]
		if len(f)==0:
			continue
		if f not in directories:
			print "No directory for "+file+"found!!"
		else:
			try :
				os.rename(file,f+'/'+file)
			except OSError:
				pass
if '-h' in sys.argv or '-help' in sys.argv :
	print '-h : help'
	print '-help : help'
	print '-m : mode for grun'
	print '-all : for running all .g16 or.inp in folder'
else: 
	if '-s' in sys.argv:
		st=sys.argv[sys.argv.index('-t')+1]
		schedule(st)
	else:
		gauss_claculations([])

#time.sleep(5)
#organize()



