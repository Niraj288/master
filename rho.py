import xlwt
import xlrd
from xlutils.copy import copy as xl_copy

filename=raw_input('Enter path : ')
file =open(filename,'r')

ref=0
d={}
n=0
for line in file:
	if 'Type' in line and 'BCP' in line:
		d[(n,'bn')]=line.strip().split()[4:]
		ref=1
	if ref==1:
		#print line
		if 'Rho' in line.strip().split():
			d[(n,'rho')]=line.strip().split()[2]
		if 'DelSqRho' in line:
			d[n,'DelSqRho']=line.strip().split()[2]
		if 'Bond Ellipticity' in line:
			d[n,'Bond Ellipticity']=line.strip().split()[3]
		if 'V' in line.strip().split():
			d[n,'V']=line.strip().split()[2]
		if 'G' in line.strip().split():
			d[n,'G']=line.strip().split()[2]
		if 'K' in line.strip().split():
			d[n,'K']=line.strip().split()[2]
		if 'L' in line.strip().split():
			d[n,'L']=line.strip().split()[2]
			ref=0
			n+=1


rb=xlrd.open_workbook('Data.xls', formatting_info=True)
#workbook = xlwt.Workbook('Data.xls')
workbook = xl_copy(rb)
name=filename.split('/')[-1]
sheet = workbook.add_sheet(name)

sheet.write(0,0,'Bond')
sheet.write(0,1,'Rho')
sheet.write(0,2,'DelSqRho')
sheet.write(0,3,'Bond Ellipticity')
sheet.write(0,4,'V')
sheet.write(0,5,'G')
sheet.write(0,6,'K')
sheet.write(0,7,'L')

index=1
max_d=-99999
for i in d:
	a,b=i
	if a>max_d:
		max_d=a
print max_d


for i in range (max_d+1):
	sheet.write(index,0,d[(index-1,'bn')])
	sheet.write(index,1,d[(index-1,'rho')])
	sheet.write(index,2,d[(index-1,'DelSqRho')])
	sheet.write(index,3,d[(index-1,'Bond Ellipticity')])
	sheet.write(index,4,d[(index-1,'V')])
	sheet.write(index,5,d[(index-1,'G')])
	sheet.write(index,6,d[(index-1,'K')])
	sheet.write(index,7,d[(index-1,'L')])
	index+=1


workbook.save('Data.xls')







