import os
import xlwt

def get_energy(file):
	f_=open(file,'r')
	f=f_.readlines()
	f_.close()
	for i in f:
		if 'FINAL SINGLE POINT ENERGY' in i:
			print i
			return float(i.split()[4])
	return 0.0

workbook = xlwt.Workbook()
sheet = workbook.add_sheet('Demo')
sheet.write(0,0,'Name')
sheet.write(0,1,'Dimer Energy')
sheet.write(0,2,'Monomer 1')
sheet.write(0,3,'Monomer 2')
sheet.write(0,4,'Delta E (Kcal/mol)')
ref=1
for i in os.listdir('.'):
	if i[-9:]=='.orca.out':
		energy=get_energy(i)
		li=i.split('.')[0].split('_')
		m1,m2=li[0],'_'.join(li[1:])
		energy_m1=get_energy('./Monomers/'+m1+'.orca.out')
		energy_m2=get_energy('./Monomers/'+m2+'.orca.out')
		deltaE=energy-(energy_m1+energy_m2)
		sheet.write(ref,0,i)
		sheet.write(ref,1,energy)
		sheet.write(ref,2,energy_m1)
		sheet.write(ref,3,energy_m2)
		sheet.write(ref,4,deltaE*627.51)
		ref+=1
workbook.save('Orca.xls')

