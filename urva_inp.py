import os
import sys

def get_s2(path):
	f=open(path,'r')
	lines=f.readlines()
	f.close()
	index,ref=-1,0
	st=''
	for line in lines:
		if 'Optimized Parameters' in line:
			index=ref
		ref+=1
	ref=0
	for i in range (index,len(lines)):
		if ref==2:
			if len(lines[i].strip().split())<3:
				pass
			else:
				para=lines[i].strip().split()[2]
				params=para[2:len(para)-1].split(',')
				st+='std    '+'    '.join(params)+'    :    '+'-'.join(para.split(','))+'\n' 
		if '-----------' in lines[i]:
			ref+=1
		if ref==3:
			break
	print st
	return st 



def string(path):
	s1="""#
# This is an input file for pURVA
#

# Control keywords:
@DATAFILETYPE = old   # old, new, xyz
@PARM = All # read in user defined parameters (No, GeomOnly, All)
@VIBRATION = off
@DIRCURV = on # calculate the reaction path direction and curvature 
@AVAM = off # calculate the adiabatic mode coupling coefficient   
@CURVCPL = off # calculate the curvature coupling coefficient 
@CORIOLIS = off # calculate the coriolis coupling coefficient 

@ENERGY = on # 0 - no ouput; 1 - energy only; 2 - energy with derivatives 
@ADIABFC = off





@DATAFILEPATH = "IRC.browse"

@BASEPATH = "./"




#Title
TITLE
   HCN
END TITLE

# Parameter list:
PARAMETER
"""
	s2=get_s2(path)
	s3="""END PARAMETER


offCURVCOR
Ln = 25
Rn = 25
END CURVCOR

offAUTOSMTH
StepSize = 0.003
Ln = 3
Rn = 3
d2ythresh = 2.4
END AUTOSMTH

offRMSPK
CutHigh = 20.0
Percentage = 0.85
GradRatio = 1.2
END RMSPK



# DMO section 
DMO
###Thresh = 1e-17
Sthresh = 0.980
Slowest = 0.890
Np = 5
NMax = 90
Cut = 0
CutA =  -0.199998264595
Skip = 1
SkipA = 1
SkipB = 1
END DMO


# CORIOLIS section
CORIOLIS

END CORIOLIS

"""
	g=open(''.join(path.split('/')[-1].split('.')[:-1])+'.inp','w')
	g.write(s1+s2+s3)
	g.close()
	return

if __name__=='__main__':
	string(sys.argv[1])





