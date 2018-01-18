import os
import subprocess
def func(path):
	
	f1=open(path,'r')
	f2=subprocess.check_output(['gcartesian', path+'.out'])
	
	file1=f1.readlines()
	file2=f2.split('\n')
	f1.close()
	
	g=open(path,'w')
	s0='%chk='+path.split('.')[0]
	s1="""
%nprocshared=8
%mem=30GB
#p counterpoise=2  TPSSTPSS empiricaldispersion=gd3bj int=ultrafine genecp

Title

0 1 0 1 0 1
"""
	g.write(s0)
	g.write(s1)
	for i in range (1,len(file2)):
		if len(file2[i])==0:
			continue
		li=file2[i].split()
		if '53'==li[0]:
			li[0]=li[0]+'(fragment=2)'
		else:
			li[0]=li[0]+'(fragment=1)'
		li='  '.join(li)
		file2[i]=li 
	s3='\n'.join(file2[1:])
	g.write(s3+'\n')
	ref=0
	for line in file1:
		if len(line.strip().split())==0 and ref==0:
			continue
		if ref>=4:
			g.write(line)
		if len(line.strip().split())==0 and ref>0:
			ref+=1
		if len(line.strip().split())>0 and '#p'== line.strip().split()[0]:
			ref=1	

	g.close()

def job(path):
	if 'g16'==path.split('.')[-1][-4:]:
		print path
		func(path)


