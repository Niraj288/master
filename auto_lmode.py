import os
import sys
import math
import subprocess
import itertools

def distance(a,b):
	a=map(float,a)
	b=map(float,b)
	res=0
	for i in range (len(a)):
		res+=(a[i]-b[i])**2
	return math.sqrt(res)

def get_log(filename):
    path=filename+'.g16.out'
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    index,ref=-1,0
    st=''
    at_ref = 0
    d = {}
    for line in lines:
        if at_ref == 1 and len(line.strip().split()) == 0:
            at_ref = 100
        if at_ref == 1:
            d[len(d)+1] = line.strip().split()[0]
        if 'Charge =' in line and at_ref == 0:
            at_ref = 1
        if 'Optimized Parameters' in line or 'Initial Parameters' in line:
            index=ref
        ref+=1
    ref=0
    
    r,a=0,0
    lm_st=''

    def get_parm_name(parm):
        p = map(int, parm)
        p = map(abs, p)

        p = map(lambda x: d[x]+str(x), p)
        print p 
        return ' : '+'-'.join(p)
    
    for i in range (index,len(lines)):
        if ref==2:
            if len(lines[i].strip().split())<3:
                pass
            else:
                para=lines[i].strip().split()[2]
                params=para[2:len(para)-1].split(',')
                #print lines[i]
                if int(params[-1])<0:
                    lm_st+=' '.join(params[:3])+' '+params[4]+get_parm_name(params[:3]+[params[4]])+'\n'
                    #print ' '.join(params),'is modified to : ',' '.join(params[:3])+' '+params[4]
                else:
                    lm_st+=' '.join(params)+get_parm_name(params)+'\n'
                if para[0]=='R':
                    r+=1
                elif para[0]=='A':
                    a+=1
                else:
                    a+=1
        if '-----------' in lines[i]:
            ref+=1
        if ref==3:
            break
    return lm_st

def get_lm_st(s):
    bl=list(itertools.combinations(s,2))
    a=list(itertools.combinations(s,3))
    #a=list(itertools.permutations(s,3))
    d=list(itertools.combinations(s,4))
    #d=list(itertools.permutations(s,4))
    st=''
    for item in bl+a+d:
        if len(item)==len(set(item)):
            st+=' '.join(item)+'\n'
    return st

def from_zmat(file):
    os.system('babel '+file+' '+file[:-5]+'.fh')
    f = open(file[:-5]+'.fh', 'r')
    lines = f.readlines()
    f.close()

    ids = int(lines[1].strip())
    count = 1
    st = ''
    for line in lines[2:]:
        k = line.strip().split()
        if len(k) < 3:
            pass
        elif len(k) < 5:
            # bond
            st += str(count) + ' ' + k[1] + '\n'
        elif len(k) < 7:
            # bond and angle
            st += str(count) + ' ' + k[1] + '\n'
            st += str(count) + ' ' + k[1] + ' '+k[3] + '\n'
        else:
            # bond, angle and dihedral
            st += str(count) + ' ' + k[1] + '\n'
            st += str(count) + ' ' + k[1] + ' '+k[3] + '\n'
            st += str(count) + ' ' + k[1] + ' '+k[3] + ' '+ k[5] + '\n'
        count += 1
    #print st
    return st



def get_cord__fchk(path):
    global atom1,atom2,d_min,d_max
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    s,cord=[],[]
    ref=0
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
            #print line.strip().split()
            lis=[i for i in line.strip().split()]
            s+=lis 
        if 'Current cartesian coordinates' in line:
            ref+=1
        if ref==4 or ref==5:
            ref+=1
        if ref==6:
            try:
                float(line.strip().split()[0])
            except ValueError:
                break
            lis=[str((float(i.split('E')[0])*(10**int(i.split('E')[1])))*0.529177) for i in line.strip().split()]
            cord+=lis
        
    s=s[:len(cord)/3]
    for i in range (len(s)):
        s[i]=str(int(float(s[i])))
    print s

    filename=path.split('/')[-1].split('.')[0]
    f=open(path.split('/')[-1].split('.')[0]+'.alm1','w')
    s1="""
$contrl
   qcprog="gaussian"
   iprint=2
   isymm = 1
   ifmatlab=.true.
   iredun = 0
 $end

$qcdata
 """
    s2='fchk="'+filename+'.fchk"\n'
    s3="""$end

$LocMod $End
"""
    lm_st=get_log(filename)
    #s4 = from_zmat(filename+'.fchk')
    s4=lm_st # to get parameters directly from gaussian output
    #s4=get_lm_st(map(str,range(1,len(s)+1))) # get permutations of all atoms

    f.write(s1+s2+s3+s4)
    f.close()
    os.system("lmode -b "+'< '+filename+'.alm1' +' >'+' '+filename+'.out1')

    

path=sys.argv[1]
dir='/'.join(path.split('/')[:-1])
if dir=='':
    pass
else:
    os.chdir(dir)
get_cord__fchk(sys.argv[1])
os.system('mv job.m '+sys.argv[1].split('.')[0]+'.m')


















   
