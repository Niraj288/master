import os 
import module
import xlwt
import subprocess


def func(path):
	print path
	global ref1,ref2
	data = subprocess.check_output('BSSE2DFT '+path, shell=True)
	data=data.split('==========================')
	if len(data[0].strip().split('\n'))<5:
		print 1234
		return
	if ('bsse' in path or 'BSSE'in path ) and 'Monomers' in path:
		if 'counterpoise=2' in ''.join(open(path,'r').readlines()[:1000]):
                        return
		energy=float(data[0].split('\n')[2].split()[4])
		conv=data[1].split('\n')[2].split()[-1]
		sheet.write(ref2,10,path.strip().split('/')[-1].split('.')[0])
		sheet.write(ref2,11,energy)
		sheet.write(ref2,12,conv)
		sheet.write(ref2,13,path)
		ref2+=1
	else:
		if 'counterpoise=2' not in ''.join(open(path,'r').readlines()[:1000]):
                	return
		rf=1
		sheet.write(ref1,0,path.strip().split('/')[-1].split('.')[0])
		lis=[]
		for i in data[0].split('\n'):
			#print i
			if len(i.strip().split())==0 or 'SCF ' not in i:
				continue
			sheet.write(ref1,rf,float(i.strip().split()[4]))
			lis.append(float(i.strip().split()[4]))
			rf+=1
		sheet.write(ref1,6,(lis[0]-(lis[3]+lis[4]))*627.51)
		sheet.write(ref1,7,(lis[0]-(lis[1]+lis[2]))*627.51)
		sheet.write(ref1,8,path)
		ref1+=1


path=raw_input("Enter path : ") or '.'

wb=xlwt.Workbook() 
sheet = wb.add_sheet(raw_input('Enter Sheet name : ') or 'bsse')
ref1=1
sheet.write(0,0,'Name')
sheet.write(0,1,'AB(AB)')
sheet.write(0,2,'A(AB)')
sheet.write(0,3,'B(AB)')
sheet.write(0,4,'A(A)')
sheet.write(0,5,'B(B)')
sheet.write(0,6,'Eint-raw\n(kcal/mol)')
sheet.write(0,7,'Eint-CP\n(kcal/mol)')
sheet.write(0,8,'comment')
sheet.write(0,10,'Name')
sheet.write(0,11,'Energy')
sheet.write(0,12,'Convergence')
sheet.write(0,13,'Comment')
ref2=1
module.search_deep(path,func,['.g16.out','.g09.out'])

wb.save('bsse_raw.xls')

