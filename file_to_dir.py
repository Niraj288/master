import os 
import module

def func(path):
	global dire
	if path==dire:
		return 
	os.system("cp "+path+' '+dire)


path=raw_input("Enter path : ")
args='['+raw_input("Enter arguments ('.g16','.g09', ...) : ")+']'
dire=path+'/'+raw_input("Enter directory name : ")+'/'
os.chdir(path)
os.system("mkdir "+dire)
module.search_deep(path,func,eval(args))