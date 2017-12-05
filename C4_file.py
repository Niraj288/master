import os
import sys
import module
def filter(lines):
	lines=lines[5:]
	lines2=[]
	for i in range(len(lines)):
		if len(lines[i].strip().split())==0:
                        break
		if 'X' in lines[i]:
			lines2.append(' '.join(lines[i].split(',')))
		elif 'Variables' in lines[i]:
			lines2.append('\n')
		else:
			li=lines[i].strip().split(',')
			for i in range (len(li)):
				if (i+1)%2!=0 and i+1!=1 and i+1<8:
					li[i]=li[i]+'*'
			lines2.append(' '.join(li)+'\n')
	return ''.join(lines2)

def make_C4(path):
    module.make_xyz(path)
    file=path.split('.')[0]+'.xyz'
    os.system('module load gaussian/g16a')
    filename=path.split('.')[0]+'.zmat'
    os.system('newzmat -ixyz -Ozmat -rebuildzmat '+file+' '+filename)
    f=open(filename,'r')
    lines=f.readlines()
    f.close()

    s1=file.split('.')[0]+".cfour : CCSD/aug-cc-pVTZ -- Comment line--\n\n"
    s2=filter(lines)
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

if __name__ == "__main__":
	make_C4(sys.argv[1])



