import os
import sys
import multiprocessing
import check_cal
import datetime
import time
import bsse_input
import threading
import nbo_input
import file_to_dir

def get_xyz(path):
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    lis=''
    ref=0
    sym={'15':'P','14':'Si','1':'H','7':'N','8':'O','6':'C','53':'I','36':'Kr','9':'F'}
    for line in lines:
        if '#p ' in line:
            ref+=1
        if len(line.strip().split())==0:
            ref+=1
        if ref==6:
            break
        if ref==3 or ref==4:
            ref+=1
        if ref==5:
            list1=line.strip().split()
            try:
                if len(line.strip().split())==2:
                    pass
                else:
                    int(list1[0])
                    lis+=sym[list1[0]]+'  '+' '.join(list1[1:])+'\n'
            except ValueError:
                lis+=line
    return lis

def threa(xyz,name,monomer,mode):
	if monomer=='A':
		f=open('/Users/47510753/Documents/master/A.txt','r')
		lines=f.readlines()
		f.close()
	else:
		lines=[]
	lines=[]
	os.chdir('./BSSE/Monomers/')
	g=open(name+'.g16','w')
	g.write('%chk='+name+'.chk\n')
	g.write("""%nprocshared=36
%mem=50GB
#p wB97XD/cc-pvTZ opt=(tight,MaxCycles=100) symmetry=loose freq int(ultrafine)

Title Card Required

0 1
""")
	g.write(xyz+'\n\n')
	g.write(''.join(lines)+'\n')
	g.write(name+'.wfx\n\n')
	g.close()
	os.system('grun '+name+'.g16'+' '+mode)

def make_bsse(name):
	print 'Making bsse for',name
	os.system('cp '+name+' BSSE/')
	os.system('cp '+name+'.out BSSE/')
	os.chdir('./BSSE')
	bsse_input.job(name)
	os.system('grun '+name+' '+mode)

def make_nbo(name):
	os.system('cp '+name+' NBO/')
	os.system('cp '+name+'.out NBO/')
	os.chdir('./NBO')
	print 'Making NBO for',name
	nbo_input.job(name)
	os.system('grun '+name+' '+mode)


def bsse_nbo(files):
	it=0
	d={}
	while it!=1:
		for i in files:
			if i not in d:
				d[i]=1
			if d[i] and check(i+'.out'):
				multiprocessing.Process(target=make_bsse,args=(i,)).start()
				multiprocessing.Process(target=make_nbo,args=(i,)).start()
				d[i]=0
		res=0
		for i in d:
			res+=d[i]
		if res==0:
			it=1
		else:
			time.sleep(5)
	return

def wfx():
	os.chdir('./wfx_files')
	os.system('python /users/nirajv/master/aimall.py')
	os.system('python /users/nirajv/master/Mwfn.py')


def check(path):
	print 'Checking status for',path,'on',datetime.datetime.now(),'\n'
	c=check_cal.func(path)
	if c==1:
		print 'Checking Done !\n'
	return c 

def submitJob(path):
	files=[]
	mode='gpu1'
	for i in os.listdir(path):
		if i[-4:]=='.g16':
			files.append(i)
			#os.system('grun '+i+' '+mode)
	os.system('mkdir BSSE')
	os.system('mkdir BSSE/Monomers')
	os.system('mkdir NBO')
	processes=[]
	for i in files:
		print i
		a,b=i[:-4].split('-')
		xyz=get_xyz(i).split('\n')
		p=multiprocessing.Process(target=threa, args=('\n'.join(xyz[:-6]),a,'A',mode,))
		processes.append(p)
		p=multiprocessing.Process(target=threa, args=('\n'.join(xyz[-6:]),b,'B',mode,))
		processes.append(p)
	for p in processes:
	        p.start()

	#threading.Thread(target=bsse_nbo, args=(files,)).start()
	bsse_nbo(files)

	file_to_dir.job('.',['wfx'],'wfx_files')
	multiprocessing.Process(target=wfx,args=()).start()

submitJob(sys.argv[1])






















