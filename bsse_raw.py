import os 
import module
import xlwt
import subprocess


def func(path):
	global ref
	data = subprocess.check_output('BSSE2DFT '+path, shell=True)
	data=data.split('==========================')
	if len(data[0].strip().split())==0:
		return
	if 'Monomers' in path:
		energy=data[0].split('\n')[0][4]
		conv=data[1].split('\n')[1].split()[-1]
		sheet.write(ref2,10,path)
		sheet.write(ref2,11,energy)
		sheet.write(ref2,12,conv)
		ref2+=1
	else:
		rf=1
		sheet.write(ref1,0,path)
		lis=[]
		for i in data[0]:
			if len(i.strip().split())==0:
				continue
			sheet.write(ref1,rf,i.strip().split()[4])
			lis.append(float(i.strip().split()[4]))
			rf+=1
		sheet.write(ref1,6,lis[0]-(lis[3]+lis[4]))
		sheet.write(ref1,7,lis[0]-(lis[1]+lis[2]))
		ref1+=1


path=raw_input("Enter path : ")

wb=xlwt.Workbook() 
sheet = wb.add_sheet(raw_input('Enter Sheet name : ') or 'pincer complexes')
ref1=1
sheet.write(0,0,'Name')
sheet.write(0,1,'AB(AB)')
sheet.write(0,2,'A(AB)')
sheet.write(0,3,'B(AB)')
sheet.write(0,4,'A(A)')
sheet.write(0,5,'B(B)')
sheet.write(0,6,'Eint-raw\n(kcal/mol)')
sheet.write(0,7,'Eint-CP\n(kcal/mol)')

sheet.write(0,10,'Name')
sheet.write(0,11,'Energy')
sheet.write(0,12,'Convergence')
ref2=1
module.search_deep(path,func,['.g16.out'])

wb.save('bsse_raw.xls')

