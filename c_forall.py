import os 
import module
import sys

def func(path):
	global args
	for arg in args:
		if '-f' in arg:
			s=arg.replace('-f',path)
			os.system(s)
			sys.stdout.flush()
		else:
			os.system(arg)
			sys.stdout.flush()

path=raw_input("Enter path : ")
lis='['+raw_input("Enter keyword for filename : ")+']'
args=eval('['+raw_input("Enter comma separated Commands : ")+']')
module.search_deep(path,func,eval(lis))
