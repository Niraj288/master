#!/usr/bin/python
import os
path=raw_input('Enter path :')
ext=raw_input('Any particular extension : ')
ext_ref=raw_input('Do you want to search only extensions (1/0) : ')
s=ext
if ext_ref!='1':
	s=raw_input('What do you want to search(including extension) : ')
k_ref=0
def search_deep(n_path):
	global k_ref
	l=len(ext)
	try:
		for i in os.listdir(n_path):
			if os.path.isdir(n_path+'/'+i):
				if i==s:
				    print ':) found a directory :'+n_path+'/'+i
				search_deep(n_path+'/'+i)
			else:
				#print i,'kdjfksjd',s[0]
				if ext_ref=='1':
					if ext==i[-l:] and l>1:
						k_ref+=1
						print '-> found extention : '+n_path+'/'+i
				elif i==s:
					k_ref+=1
					print '-> found a file : '+n_path+'/'+i 
				elif i==s+ext and l>1:
					k_ref+=1
					print '-> found a file : '+n_path+'/'+i 
				#search(i,n_path)
	except OSError:
		raise Exception('Not a directory ')
search_deep(path)

if k_ref==0:
	print "\n:( no files found in the specified path :'( "


