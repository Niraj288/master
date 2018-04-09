import module
import os

def func(path):
	global d
	f=open(path,'r')
	lines=f.readlines()[-10:]
	f.close()
	k=[]
	for line in lines:
		if 'Job cpu time' in line:
			k=line.strip().split()
			d['d']+=float(k[3])
			d['h']+=float(k[5])
			d['m']+=float(k[7])
			d['s']+=float(k[9])
	print path,k[3:10]
	return

d={'d':0,'h':0,'m':0,'s':0}
path='/'.join(os.getcwd().split('/')[:3])+'/test'
print 'Searching in',path
module.search_deep(path,func,['.out'])

time=d['d']*24*60*60+d['h']*60*60+d['m']*60+d['s']

day = time // (24 * 3600)
time = time % (24 * 3600)
hour = time // 3600
time %= 3600
minutes = time // 60
time %= 60
seconds = time

print 'Days : ',day
print 'Hours : ',hour
print 'Minutes : ',minutes
print 'Seconds : ',seconds

