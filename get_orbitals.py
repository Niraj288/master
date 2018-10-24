import sys

def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush() 

def get_atoms(lines,index):
	dic={}
	for i in range (index,len(lines)):
		line=lines[i]
		if len(line.strip().split())==0:
			break
		a=line.strip().split()[0]
		if a not in dic:
			dic[a]=1
		else:
			dic[a]+=1
	return dic

def orbitals(file,d):
	f=open(file,'r')
	lines=f.readlines()
	f.close()
	ref_p=0
	#print len(lines)
	if len(lines)>10000:
		ref_p=1
	ref=0
	count=0
	sym_table,num_table=[],[]
	name,id=None,None
	total=len(lines)
	string=''
	for i in range (len(lines)):
		if 'Symbolic Z-matrix:' in lines[i]:
			atoms=get_atoms(lines,i+2)
			#print atoms
		if ref_p:
			progress(i, total, status='Progress ')
		if 'Condensed to atoms (all electrons):' in lines:
			break
		if i==ref:
			continue
		if 'Eigenvalues --' in lines[i]:
			sym_table=lines[i-1].strip().split()
			num_table=lines[i-2].strip().split()
			string+=lines[i-2]
			string+=lines[i-1]
		if 'Molecular Orbital Coefficients:' in lines[i]:
			ref=i+1
			count=1
      
		if count:
			for j in atoms:
				if j in lines[i].strip().split():
					name=j 
					id=lines[i].strip().split()[1]
					#print name,id
			if name in d and id==d[name]['id']:
				lis=lines[i].strip().split()
				for j in d[name]:
					if len(j.split())>1:
						check_li=lines[i]
					else:
						check_li=lis 
					if j in check_li:
						#count+=1
						#sys.stdout.write("\rPoints found : %i" % count)
						#sys.stdout.flush()
						string+=lines[i]#+'    '+name+' '+id
						d[name][j].append([lis[0]]+lis[-len(sym_table):]+num_table+sym_table)
	print string
	'''
	print '\n'*2
	#print d
	lale,la='',5
	for i in d:
		for j in d[i]:
			if j!='id':

				for k in d[i][j]:
					if k[-la:]==lale[-la:]:
						if len(k)<16:
							l=(len(k)-1)/3
							print k[0],i,d[i]['id'],j,' '*15+"        ".join(k[1:l+1]) 
						else:
							print k[0],i,d[i]['id'],j,'      '+"{:>15} {:>15} {:>15} {:>15} {:>15}".format(*k[1:6]) 
					else:
						if len(k)<16:
							l=(len(k)-1)/3
							print ' '*30+"             ".join(k[1+l:1+2*l])
							print ' '*30+"             ".join(k[1+2*l:])
							print k[0],i,d[i]['id'],j,' '*15+"        ".join(k[1:l+1]) 
						else:
							print ' '*15+"{:>15} {:>15} {:>15} {:>15} {:>15}".format(*k[6:11]) 
							print ' '*15+"{:>15} {:>15} {:>15} {:>15} {:>15}".format(*k[11:])
							print k[0],i,d[i]['id'],j,'      '+"{:>15} {:>15} {:>15} {:>15} {:>15}".format(*k[1:6]) 
					#print k[0],i,d[i]['id'],j,' '.join(k)
					la=(len(k)-1)/3
					lale=k[-la:]
	'''
					


def job():
	d={}
	atoms=raw_input('Enter atom and orbitals (Eg. C 1 2PX 2PY.. ) :')
	while len(atoms)!=0:
		lis=atoms.split()
		d[lis[0]]={'id':lis[1]}
		for j in lis[2:]:
			s=j.replace('_',' ')
			d[lis[0]][s]=[]

		atoms=raw_input('Enter another atom and orbitals (leave blank if not) :')
	print d
	st=orbitals(sys.argv[1],d)

job()
