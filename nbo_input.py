import os
import subprocess

def func(path):
        filename=path.split('/')[-1].split('.')[0]
        s1="#p pop=(nbo6,savenbos) symmetry=loose TPSSTPSS empiricaldispersion=gd3bj int=ultrafine genecp\n"
        s2=subprocess.check_output('gcartesian '+path+'.out', shell=True)+'\n'
        file=open(path,'r')
        lines=file.readlines()
        file.close()
        ref=0
        inde=None
        for i in range (len(lines)):
                if len(lines[i].split())>0 and '.wfx'==lines[i].strip().split()[0][-4:]:
                        lines[i]=''
                if 'procshared' in lines[i]:
                        lines[i]='%nprocshared=4\n'
                if '%mem' in lines[i]:
                        lines[i]='%mem=10GB\n'
                if '#p' in lines[i]:
                        lines[i]=s1
                        ref=1
                if ref>0 and len(lines[i].split())==0:
                        if ref==5:
                                lines[i]=''
                                inde=i+1
                        ref+=1
                if ref==3 or ref==4:
                        ref+=1
                if ref==5:
                        lines[i]=''
        g=open(path[:-4]+'.g09','w')
        lines.insert(inde,s2)
        for line in lines:
                g.write(line)
        g.close()

def job(path):
        func(path)

#func('I2_NCH32_Pt.g16')
if __name__=='__main__':
        for i in os.listdir('.'):
                if i[-4:] in ['.g16','.g09']:
                        print i
                        func(i)