#!/usr/bin/python
import os
import module
import sys
import multiprocessing
import copy

def func(path):
        global args1,path_lis
        filename=path.split('/')[-1].split('.')[0]
        ref = 0
        if '-d' in args1:
                path_lis.append(path)
                return
        #print temp_args1
        #print os.getcwd()
        for arg in args1:
                if '-p' in arg or '-f' in arg:
                        s=arg.replace('-f',filename)
                      
                        s=s.replace('-p',path)
                        print 'Performing :',s,'...'
                        os.system(s)
                        sys.stdout.flush()
                else:
                        os.system(arg)
                        sys.stdout.flush()
        

path=raw_input("Enter path : ")
o_path = os.getcwd()
path_lis = []
lis=raw_input("Enter keyword for filename : ")
lis=lis.replace(',','","')
lis=eval('["'+lis+'"]')
print '-f is for filename'
print '-p is for path'
print '-d to work in the same directory'
args1=raw_input("Enter comma separated Commands : ")
args1=args1.replace(',','","')
args1=eval('["'+args1+'"]')
module.search_deep(path,func,lis)

print path_lis

if len(path_lis)>0:
        args1.remove('-d')
        for path in path_lis:
                #print 'Back on :',os.getcwd()
                filename=path.split('/')[-1].split('.')[0]
                file_path = '/'.join(path.split('/')[:-1])
                file = path.split('/')[-1]
                os.chdir(file_path)
                #print 'working on :',os.getcwd()
                for arg in args1:
                        if '-p' in arg or '-f' in arg:
                                s=arg.replace('-f',filename)

                                s=s.replace('-p',file)
                                print 'Performing :',s,'...'
                                os.system(s)
                                sys.stdout.flush()
                        else:
                                os.system(arg)
                                sys.stdout.flush()
                os.chdir(o_path)












