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
	li=filename.split('_')
	m,m_=li[0],'_'.join(li[1:])
	mg=open('monomers'+'/'+m+'.txt')
	g=''
	for line in mg:
		if 'SCF Done' in line:
			g=float(line.strip().split()[4])

	fi=open(path,'r')
	li=[]
	for line in fi:
		if 'SCF Done' in line:
			li.append(line.strip().split()[4])
	li=map(float,li)

	#print li	
	f=''
	mf=open('monomers/'+m_+'.txt')
        for line in mf:
                if 'SCF Done' in line:
                        f=float(line.strip().split()[4])
	a,b,c,d,e=li
	Eint=a-(b+c)
	ErelaxA=d-g
	ErelaxB=e-f
	ErelaxT=ErelaxA+ErelaxB
	deltaE=Eint+ErelaxT
	print deltaE
	sheet.write(ref,0,name)
	sheet.write(ref,1,ErelaxA)
	sheet.write(ref,2,ErelaxB)
	sheet.write(ref,3,ErelaxT)
	sheet.write(ref,4,Eint)
	sheet.write(ref,5,deltaE)
	ref+=1

path=raw_input("Enter path : ")

wb=xlwt.Workbook() 
sheet = wb.add_sheet(raw_input('Enter Sheet name : ') or 'pincer complexes')
ref=1
sheet.write(0,0,'Name')
sheet.write(0,1,'Erelax_A')
sheet.write(0,2,'Erelax_B')
sheet.write(0,3,'Erelax_T')
sheet.write(0,4,'Eint')
sheet.write(0,5,'deltaE')
module.search_deep(path,func,['.txt'])

wb.save('bsse.xls')


