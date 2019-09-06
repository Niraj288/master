import sys
def sub_s(lis):
        sub_s,ref='',0
        for i in lis:
                if ref%5==0:
                        sub_s+='\n    '+i+'    '
                else:
                        sub_s+=i+'    '
                ref+=1
        return sub_s


def xyz_(lines):
	atd = {'Ru': 44, 'Re': 75, 'Rf': 104, 'Rg': 111, 'Ra': 88, 'Rb': 37, 'Rn': 86, 'Rh': 45, 'Be': 4, 'Ba': 56, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'Og': 118, 'H': 1, 'P': 15, 'Os': 76, 'Es': 99, 'Hg': 80, 'Ge': 32, 'Gd': 64, 'Ga': 31, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'C': 6, 'Pb': 82, 'Pa': 91, 'Pd': 46, 'Cd': 48, 'Po': 84, 'Pm': 61, 'Hs': 108, 'Ho': 67, 'Hf': 72, 'K': 19, 'He': 2, 'Md': 101, 'Mg': 12, 'Mc': 115, 'Mo': 42, 'Mn': 25, 'O': 8, 'Mt': 109, 'S': 16, 'W': 74, 'Zn': 30, 'Eu': 63, 'Zr': 40, 'Er': 68, 'Nh': 113, 'Ni': 28, 'No': 102, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10, 'Np': 93, 'Fr': 87, 'Fe': 26, 'Fl': 114, 'Fm': 100, 'B': 5, 'F': 9, 'Sr': 38, 'N': 7, 'Kr': 36, 'Si': 14, 'Sn': 50, 'Sm': 62, 'V': 23, 'Sc': 21, 'Sb': 51, 'Sg': 106, 'Se': 34, 'Co': 27, 'Cn': 112, 'Cm': 96, 'Cl': 17, 'Ca': 20, 'Cf': 98, 'Ce': 58, 'Xe': 54, 'Lu': 71, 'Cs': 55, 'Cr': 24, 'Cu': 29, 'La': 57, 'Ts': 117, 'Li': 3, 'Lv': 116, 'Tl': 81, 'Tm': 69, 'Lr': 103, 'Th': 90, 'Ti': 22, 'Te': 52, 'Tb': 65, 'Tc': 43, 'Ta': 73, 'Yb': 70, 'Db': 105, 'Dy': 66, 'Ds': 110, 'I': 53, 'U': 92, 'Y': 39, 'Ac': 89, 'Ag': 47, 'Ir': 77, 'Am': 95, 'Al': 13, 'As': 33, 'Ar': 18, 'Au': 79, 'At': 85, 'In': 49}
	count, c, m, xyz , amass, za = 0, None, None, [], [], []
	for line in lines:
		count += 1
		if 'Geometry (in Angstrom)' in line:
			k = line.strip().split()
			c, m = k[5][:-1], k[8][:-1]
		if 'Center              X                  Y                   Z ' in line:
			for j in range (count+1, len(lines)):
				if len(lines[j].strip().split())==0:
					break
				kl = lines[j].strip().split()
				za.append(str(atd[kl[0]]))
				amass.append(kl[-1])
				xyz += kl[1:4]
			break
	for i in range (max(0, len(lines)-1000), len(lines)):
		if 'Force constants in Cartesian coordinates.' in lines[i]:
			break
	hess = []
	#print (lines[i])
	for line in lines[i+2:]:
		if len(line.strip().split()) == 0:
			break
                hess += line[4:].strip().split()
	for i in range (len(hess)):
		hess[i] = hess[i].replace('[', '')
		hess[i] = hess[i].replace(']', '')
	return xyz , amass, za, hess
	



def make_dat(path):
        filename=path.split('/')[-1].split('.')[0]
        #print filename[:-4]+'_hessian'
	f = open(path, 'r')
	lines = f.readlines()
	f.close()
	xyz , amass, za, hess = xyz_(lines)
        g=open(filename+'.dat','w')
        s=''
        s+='Eventually the world will end in singularity !! \n'
        s+='NATM\n   '+str(len(za))+'\nAMASS'+sub_s(amass)
        s+='\nZA'+sub_s(za)
        s+='\nXYZ'+sub_s(xyz)
        s+='\nNOAPT'
        s+='\nFFX'+sub_s(hess)
        g.write(s)
        g.close()

	h = open(filename+'.alm', 'w')
	h.write('''
 $contrl
   qcprog="ALMODE"
   iprint=0
   isymm = 1
   ifsave=.false.
 $end

$qcdata
  fchk="''' + filename + '''.dat"
$end

$LocMod $End
''')
	h.close()



if __name__ == '__main__':
	make_dat(sys.argv[1])
