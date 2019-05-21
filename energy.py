import os
import subprocess
import xlwt

lis=[]
for i in os.listdir('.'):
	if i[-4:]=='.out':
		data = subprocess.check_output("grep 'SCF Done' "+ i +' | tail -1', shell=True)
		data2 = subprocess.check_output("grep 'Sum of electronic and thermal Enthalpies' "+ i +' | tail -1', shell=True)
		data3 = subprocess.check_output("grep 'Sum of electronic and thermal Free Energies' "+ i +' | tail -1', shell=True)

		if len(data)<5:
			pass
		else:
			if len(data2) > 5:
				lis.append([i,data.split('\n')[0].split()[4],data2.strip().split()[-1],data3.strip().split()[-1]])
			else:
				lis.append([i,data.split('\n')[0].split()[4],'-','-'])


wb=xlwt.Workbook() 
sheet = wb.add_sheet('Energy')

lis.sort()
sheet.write(0,0,'Filename')
sheet.write(0,1,'Energy (Hartree)')
sheet.write(0,2,'Enthalpy (Hartree)')
sheet.write(0,3,'Gibbs Free Energy (Hartree)')

for i in range (len(lis)):
	sheet.write(i+1,0,lis[i][0])
	sheet.write(i+1,1,lis[i][1])
	sheet.write(i+1,2,lis[i][2])
	sheet.write(i+1,3,lis[i][3])

wb.save('Energy.xls')




