#!/usr/bin/python
import xlwt
import module
#import xlrd
#from xlutils.copy import copy as xl_copy

ref=0
d={}
n=0
metals=['Co','Rh','Ni','Pd','Pt','Ir','I']
def func(path):
        global n,d,ref
        filename=path
        file =open(filename,'r')
        for line in file:
                if 'Type' in line and 'BCP' in line:
                        d[(n,'bn')]=line.strip().split()[4:]
                        d[(n,'comment')]=path
			ref=1
                if ref==1:
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
        file.close()

def check(string):
	a,b=string
	if 'Ir' in a or 'Ir' in b:
		if 'I' in a and 'I' in b:
			return True
		else:
			return False
        if 'I' in a or 'Cl' in a or 'Br' in a:
                return True
	if 'I' in b or 'Cl' in b or 'Br' in b:
		return True
        return False

path=raw_input("Enter path : ")
module.search_deep(path,func,['.mgpvig'])


#rb=xlrd.open_workbook('Data.xls', formatting_info=True)
workbook = xlwt.Workbook()
#workbook = xl_copy(rb)
name=raw_input('Enter sheet name : ')
sheet = workbook.add_sheet(name)

sheet.write(0,0,'Bond')
sheet.write(0,1,'Rho')
sheet.write(0,2,'DelSqRho')
sheet.write(0,3,'Bond Ellipticity')
sheet.write(0,4,'V')
sheet.write(0,5,'G')
sheet.write(0,6,'K')
sheet.write(0,7,'L')
sheet.write(0,8,'Comments')
index=1
max_d=-99999
for i in d:
        a,b=i
        if a>max_d:
                max_d=a
print max_d
i_ref=1
print d
for i in range (max_d+1):
        if check(d[(index-1,'bn')]):
		d[(index-1,'bn')]='-'.join(d[(index-1,'bn')])
                sheet.write(i_ref,0,d[(index-1,'bn')])
                sheet.write(i_ref,1,d[(index-1,'rho')])
                sheet.write(i_ref,2,d[(index-1,'DelSqRho')])
                sheet.write(i_ref,3,d[(index-1,'Bond Ellipticity')])
                sheet.write(i_ref,4,d[(index-1,'V')])
                sheet.write(i_ref,5,d[(index-1,'G')])
                sheet.write(i_ref,6,d[(index-1,'K')])
                sheet.write(i_ref,7,d[(index-1,'L')])
		sheet.write(i_ref,8,d[(index-1,'comment')])
        	i_ref+=1        
	index+=1


workbook.save('Data.xls')
