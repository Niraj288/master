import module
import os
from datetime import datetime, timedelta

def log_file(path):
	m_d={'Jan':'01','Feb':'02','Mar':'03','Apr':'04','May':'05','Jun':'06','Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
	f=open(path,'r')
	lines=f.readlines()
	f.close()
	then,now='',''
	for line in lines[:100]:
		if '-- Start Module:' in line:
			
			month,day,time,year=line.strip().split()[6:10]
			st = year+'-'+m_d[month]+'-'+(2-len(day))*'0'+day+'T'+time+'Z'
			then = datetime.strptime(st, '%Y-%m-%dT%H:%M:%SZ')
			break
	
	lines=lines[-20:]
	for line in lines:
		if '-- Stop Module:' in line:
			month,day,time,year=line.strip().split()[6:10]
			st = year+'-'+m_d[month]+'-'+(2-len(day))*'0'+day+'T'+time+'Z'
			now = datetime.strptime(st, '%Y-%m-%dT%H:%M:%SZ')
			break
	if then=='':
		return [0,0,0,0]
	else:
		ft = str(now-then)
		day=int(ft[0])
		h,m,s=map(float,ft.split()[-1].split(':'))
		return [day,h,m,s]

def func(path):
	if './.'==path[:3]:
		return
	global d
	f=open(path,'r')
	lines=f.readlines()[-10:]
	f.close()
	k=[]
	if path[-4:]=='.log':
		if 'MOLCAS' in ''.join(lines[:50]):
			k=log_file(path)
			d['d']+=k[0]
			d['h']+=k[1]
			d['m']+=k[2]
			d['s']+=k[3]
			print path,k
			return

	for line in lines:
		if 'Job cpu time' in line:
			k=line.strip().split()
			if k[9]=='****':
				k[9]='0.0'
			d['d']+=float(k[3])
			d['h']+=float(k[5])
			d['m']+=float(k[7])
			d['s']+=float(k[9])
			print path,k[3:10]
		elif 'CPU TIMES  *' in line:
			k=line.strip().split()
			d['s']+=float(k[3].split('.')[0])
			print path,[0,0,0,k[3]]
		
	return

d={'d':0,'h':0,'m':0,'s':0}
path=raw_input("Enter path : ") or '/'.join(os.getcwd().split('/')[:3])+'/scratch'
print 'Searching in',path
module.search_deep(path,func,['.out','.log'])

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

