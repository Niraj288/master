import os
import sys
import module

def make_C4(path):
	module.make_out(path)
	file=path.split('.')[0]+'.xyz'
	os.system('module load gaussian/g16a')
    filename=path.split('.')[0]+'.zmat'
    os.system('newzmat -ixyz -Ozmat -rebuildzmat '+file+' '+filename)
    f=open(filename,'r')
    lines=f.readlines()
    f.close()

    s1=file.split('.')[0]+".cfour : CCSD/aug-cc-pVTZ -- Comment line--\n"
    s2=''.join(lines)
    s3="""
*ACES2(CALC=CCSD(T)
ABCDTYPE=AOBASIS
CC_PROG=ECC
BASIS=AUG-PVTZ
FROZEN_CORE=ON
SCF_CONV=10
SCF_MAXCYC=300
CC_CONV=10
CC_MAXCYC=200
LINEQ_CONV=10
LINEQ_MAXCY=150
OPT_MAXCYC=150
GEO_CONV=7
MEM_UNIT=GB
MEMORY=60
CHARGE=0
MULTIPLICITY=1)
"""
	g=open(path.split('.')[0]+'.cfour','w')
	g.write(s1+s2+s3)
	g.close()

if len(sys.argv)>1:
	lis=[sys.argv[1]]
else:
	lis=os.listdir('.')

for i in lis:
	if i[-4:]=='.out':
		make_C4(i)