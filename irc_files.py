import os
import sys

name='.'.join(sys.argv[1].split('.')[:-1])

st1='%chk='+name+"""_for.chk
%NProcShared=28
%mem=132GB

#P BP86/6-31G(d,p) scf(conver=12) irc(calcall,tight,maxcycle=50,stepsize=-3,maxpoints=999,forward) guess(read) int(acc2e=12,SuperFineGrid) geom(allcheck) iop(1/40=8) iop(1/108=-1,1/109=-1) iop(1/45D=1000000) iop(1/131=12201)

"""+name+'_for.urv\n\n'

st2='%chk='+name+"""_rev.chk
%NProcShared=28
%mem=40GB

#P BP86/6-31G(d,p) scf(conver=12) irc(calcall,tight,maxcycle=50,stepsize=-3,maxpoints=999,reverse) guess(read) int(acc2e=12,SuperFineGrid) geom(allcheck) iop(1/40=8) iop(1/108=-1,1/109=-1) iop(1/45D=1000000) iop(1/131=12201)

"""+name+'_rev.urv\n\n'

os.system('cp '+sys.argv[1]+' '+name+'_for.chk')
os.system('cp '+sys.argv[1]+' '+name+'_rev.chk')

f=open(name+'_for.g09','w')
f.write(st1)
f.close()

g=open(name+'_rev.g09','w')
g.write(st2)
g.close()

st=raw_input("Enter mode : ")
if len(st)>0:
	os.system('grun '+name+'_for.g09 '+st)
	os.system('grun '+name+'_rev.g09 '+st)


