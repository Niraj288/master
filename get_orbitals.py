import sys

def orbitals(file,d):
	f=open(file,'r')
	lines=f.readlines()
	f.close()
	ref=0
	sym_table,num_table=[],[]
	name,id='',''
	for i in range (len(lines)):
		if i==ref:
			continue
		if 'Eigenvalues --' in lines[i]:
			sym_table=lines[i-1].strip().split()
			num_table=lines[i-2].strip().split()
		if 'Molecular Orbital Coefficients:' in lines[i]:
			ref=i+1

		for j in d:
			if j in lines[i].strip().split():
				name=j 
				id=lines[i].strip().split()[1]
		if name in d and id==d[name]['id']:
			lis=lines[i].strip().split()
			for j in d[name]:
				if j in lis:
					d[name][j].append(lis[-5:]+num_table+sym_table)
	for i in d:
		for j in d[i]:
			if j!='id':
				for k in d[i][j]:
					print i,d[i]['id'],j,' '.join(k)


def job():
	d={}
	atoms=raw_input('Enter atom and orbitals (Eg. C 1 2PX 2PY.. ) :')
	while len(atoms)!=0:
		lis=atoms.split()
		d[lis[0]]={'id':lis[1]}
		for j in lis[2:]:
			d[lis[0]][j]=[]

		atoms=raw_input('Enter another atom and orbitals (leave blank if not) :')
	st=orbitals(sys.argv[1],d)

job()
