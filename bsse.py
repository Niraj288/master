import os 
import module
import xlwt
#import xlrd

def func(path):
	global ref
	if 'monomers' in path:
		return
	name=path.replace('/','_')
	filename=path.split('/')[-1].split('.')[0]
	m=filename.split('_')[0]
	mf=open('monomers'+'/'+m+'.txt')
	g=''
	for line in mf:
		if 'SCF Done' in line:
			g=float(line.strip().split()[4])

	fi=open(path,'r')
	li=[]
	for line in fi:
		if 'SCF Done' in line:
			li.append(line.strip().split()[4])
	li=map(float,li)

	print li	
	sheet.write(ref,6,g)
	sheet.write(ref,0,name)
	for i in range (5):
		sheet.write(ref,i+1,li[i])
	ref+=1

path=raw_input("Enter path : ")

wb=xlwt.Workbook() 
sheet = wb.add_sheet(raw_input('Enter Sheet name : ') or 'pincer complexes')
ref=1
sheet.write(0,0,'Name')
sheet.write(0,1,'a')
sheet.write(0,2,'b')
sheet.write(0,3,'c')
sheet.write(0,4,'d')
sheet.write(0,5,'e')
sheet.write(0,6,'g')
module.search_deep(path,func,['.txt'])

wb.save('bsse.xls')


