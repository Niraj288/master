import os
def search_deep(n_path,func,args=[],ref=1):
        if 'Trash' in n_path:
		return
	try:
                for i in os.listdir(n_path):
                        if os.path.isdir(n_path+'/'+i):
                                search_deep(n_path+'/'+i,func,args,ref)
                        else:
				if ref:
                                	for j in args:
                                        	if j==i[-len(j):]:
                                                	func(n_path+'/'+i)
				else:
					for j in args:
                                                if j in i:
                                                        func(n_path+'/'+i)
        except OSError:
                 func(n_path)
                
                #raise Exception('Not a directory ')

def make_fchk(path):
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    s,cord=[],[]
    ref=0
    sym={'14':'Si','1':'H','7':'N','8':'O','6':'C','53':'I','36':'Kr','9':'F'}
    for line in lines:
        if 'Number of symbols in /Mol/' in line:
            break
        if len(line.strip().split())==0:
            continue
        if 'Nuclear charges' in line:
            ref=1
        if ref==1 or ref==2:
            ref+=1
        if ref==3:
            lis=[i.split('.')[0] for i in line.strip().split()]
            s+=lis 
        if 'Current cartesian coordinates' in line:
            ref+=1
        if ref==4 or ref==5:
            ref+=1
        if ref==6:
            lis=[str((float(i.split('E')[0])*(10**int(i.split('E')[1])))*0.529177) for i in line.strip().split()]
            cord+=lis
        
    cords=[]
    i,ref=0,0
    g=open(path[:-5]+'.xyz','w')
    while i < (len(cord)):
        g.write(sym[s[ref]]+' '+' '.join(cord[i:i+3])+'\n')
        i+=3
        ref+=1
    g.close()

def make_out(path):
    data = subprocess.check_output('gcartesian '+path, shell=True)
    g=open(path[:-4]+'.xyz','w')
    g.write(data)
    g.close()

def make_xyz(path):
    if path.split('.')[-1]=='fchk':
        make_fchk(path)
        return 2
    elif path.split('.')[-1]=='out':
        make_out(path)
        return 3
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    lis=''
    ref=0
    sym={'15':'P','14':'Si','1':'H','7':'N','8':'O','6':'C','53':'I','36':'Kr','9':'F'}
    for line in lines:
        if '#p ' in line:
            ref+=1
        if len(line.strip().split())==0:
            ref+=1
        if ref==6:
            break
        if ref==3 or ref==4:
            ref+=1
        if ref==5:
            list1=line.strip().split()
            try:
                if len(line.strip().split())==2:
                    lis+=line
                else:
                    int(list1[0])
                    lis+=sym[list1[0]]+'  '+' '.join(list1[1:])+'\n'
            except ValueError:
                lis+=line
    g=open(path[:-4]+'.xyz','w')
    g.write(lis)
    g.close()
    return 1

