import os 
import module

def func(path,args):
	os.chdir(path)
	for arg in args:
		os.system(arg)


path=raw_input("Enter path : ")
args='['+raw_input("Enter comma separated Commands : ")+']'
module.search_deep(path,func,eval(args))