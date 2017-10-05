#!/usr/bin/python
import os
try:
	import pyPdf
except ImportError:
	print 'let me install the pdf reader for ya..'
	os.system('pip install pyPdf')
	import pyPdf
path=raw_input('Enter path :')
s=raw_input('What do you want to search :')
d_ref=raw_input('Files/Directories u want to ignore :')
e_ref=raw_input('Any particular file extensions :')
d_ref=d_ref.split(',')
s=s.split()
l_s=len(s)
for i in range (l_s):
	s[i]=s[i].lower()
k_ref=0
def getPDFContent(path,i):
    content = ""
    # Load PDF into pyPDF
    pdf = pyPdf.PdfFileReader(file(path, "rb"))
    li=[]
    for page in pdf.pages:
    	li+=page.extractText().split()
    li=[x.lower() for x in li]
    path=path.split('/')
    path.remove(i)
    path='/'.join(path)
    prin(li,i,path)
def prin(lis,jk,path):
	global k_ref,l_s
	for j in range (len(lis)):
		if s[0]==lis[j]:
			if s==lis[j:j+l_s]:
				k_ref=1
				print ':) found in '+path+'/'+jk
				break

def search(i,n_path):
	#for i in os.listdir(path):
	#print n_path,'search'
	try:
		os.chdir(n_path)
		if i[-4:]=='.pdf':
			getPDFContent(n_path+'/'+i,i)
		else:
			if i in d_ref:
				pass
			else:
				k=open(i,'r')
				lis=[]
				for j in k:
					lis+=[x.lower() for x in j.split()]
				#d[i]=lis
				prin(lis,i,n_path)
	except IOError:
		print '-->>:( '+i,' is out of scope for python and so its ignored'
	
n_path=path
def search_deep(n_path):

	try:
		for i in os.listdir(n_path):
			if os.path.isdir(n_path+'/'+i):
				if i in d_ref:
					pass
				else:
				    search_deep(n_path+'/'+i)
			else:
				if len(e_ref)>0 and i[-len(e_ref):]==e_ref:
				    search(i,n_path)
				elif len(e_ref)==0:
					search(i,n_path)
	except OSError:
		raise Exception('Not a directory ')
search_deep(n_path)

if k_ref==0:
	print "not found in the specified directory :( :( "

