
import os
import xlwt

d={} 
def lmode(path):
	print path
	file = open(path,'r')
	filename=path.split('/')[-1].split('.')[0]
	global d
	filename=path 
	d[filename]=[]
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
			print line.strip().split()[5:8]
			bond,q_n,ka,wa=line.strip().split()[5:8]+[line.strip().split()[9]]
			d[filename].append(['bond,q_n,ka,wa',bond,float(q_n),float(ka),float(wa)])
		if ref_l==4:
			ref_l=0
		
	d[filename].append(['wo',''])
	d[filename].append(['path',path])




#lmode(file,path)

path=raw_input('Enter path :')

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
			sheet.write(c,5,j[1])
		else:
			sheet.write(num,1,j[1])
			sheet.write(num,2,j[2])
			sheet.write(num,3,j[3])
			sheet.write(num,4,j[4])
			num+=1
wb.save('lmode.xls')
