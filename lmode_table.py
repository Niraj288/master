#!/usr/bin/python
import os
import xlwt
import numpy as np
d={}
f={}
def lmode(path):
	print path
	file = open(path,'r')
	filename=path.split('/')[-1].split('.')[0]
	global d
	global f
	filename=path 
	d[filename]=[]
	f[filename]=[]
	ref_l=0
	ref_o,lm=0,{}
	p_ref1,p_ref2=0,0
	for line in file:
		p_ref1+=1
		if 'Program LOCALMODES' not in line and p_ref1>5 and p_ref2==0:
			return
		elif 'Program LOCALMODES' in line:
			p_ref2=1
		if 'Local mode properties:' in line:
			ref_l=1
		if '------------------' in line and ref_l>0:
			ref_l+=1
		elif ref_l==3:
			index=5
			try:
				float(line.strip().split()[6])
			except ValueError:
				index=6
			if len(line.strip().split()[index:])==6:
				st=line.strip().split()
				print [st[index]]+[st[index+1]]+[st[index+3]]+[line.strip().split()[-1]]
				bond,q_n,ka,wa=[st[index]]+[st[index+1]]+[st[index+3]]+[line.strip().split()[-1]]
                        	d[filename].append(['bond,q_n,ka,wa',bond,float(q_n),float(ka),float(wa)])
                        	f[filename].append([bond,float(q_n),float(ka),float(wa)])
			else:
				print line.strip().split()[index:index+3]+[line.strip().split()[-1]]
				bond,q_n,ka,wa=line.strip().split()[index:index+3]+[line.strip().split()[-1]]
				d[filename].append(['bond,q_n,ka,wa',bond,float(q_n),float(ka),float(wa)])
				f[filename].append([bond,float(q_n),float(ka),float(wa)])
		if ref_l==4:
			ref_l=0
		
	d[filename].append(['wo',''])
	d[filename].append(['path',path])


#lmode(file,path)

path='.' or raw_input('Enter path :')

k_ref=0
def search_deep(n_path):
	global k_ref
	try:
		for i in os.listdir(n_path):
			if os.path.isdir(n_path+'/'+i):
				search_deep(n_path+'/'+i)
			else:
				if '.out' in i:
					k_ref+=1
					lmode(n_path+'/'+i)
	except OSError:
		k_ref+=1
		lmode(n_path)
		#raise Exception('Not a directory ')
search_deep(path)
np.save('Data.npy',f)
print d

if k_ref==0:
	print "\n:( no files found in the specified path :'( "

wb=xlwt.Workbook() 

sheet = wb.add_sheet(raw_input('Enter output sheet name : ') or path.split('/')[-1][:-4])

sheet.write(0,0,'Name')
sheet.write(0,1,'Bond')
sheet.write(0,2,'Bond length')
sheet.write(0,3,'Ka')
sheet.write(0,4,'wa')
sheet.write(0,5,'Comments')
sheet.write(0,6,'Path')

num=0
for i in d:
	if len(d[i])==0:
		continue
	else:
		num+=1
	sheet.write(num,0,i.split('/')[-1].split('.')[0])
	c=num
	for j in d[i]:
		if 'wo' in j:
			continue
		elif 'path' in j:
			sheet.write(c,6,j[1])
		else:
			sheet.write(num,1,j[1])
			sheet.write(num,2,j[2])
			sheet.write(num,3,j[3])
			sheet.write(num,4,j[4])
			sheet.write(num,5,j[5])
			num+=1
wb.save('lmode.xls')
