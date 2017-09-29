#Automation
import os
import threading
import multiprocessing
import time
import sys
directories={}

def threa(command):
        os.system(command)

def gauss_claculations():
	global directories
	path,filename='',''
	if '-m' in sys.argv:
		mode=sys.argv[2]
	else:
		mode=raw_input("Enter mode : ")
	processes=[]
	refe=0
	if '-all' in sys.argv:
		refe=1
	elif raw_input('Send all .g16 or .inp files in the directory ?? (y/n) : ')=='y':
		refe=1
	for i in os.listdir(os.getcwd()):
	        if '.g16' in i or '.inp' in i:
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
if '-h' or '-help' in sys.argv:
	print '-h : help'
	print '-help : help'
	print '-m : mode for grun'
	print '-all : for running all .g16 or.inp in folder'
else: 
	gauss_claculations()
#time.sleep(5)
#organize()



