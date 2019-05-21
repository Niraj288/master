import sys

file = sys.argv[1]

f = open(file+'.job', 'w')

st = """#!/bin/bash
#SBATCH -J file
#SBATCH -o file.out
#SBATCH -N proc
#SBATCH --exclusive
#SBATCH --mem=100000
#SBATCH -p mod

export XTBHOME=/users/nirajv/xtb_exe


export OMP_NUM_THREADS 36
export MKL_NUM_THREADS 36
export OMP_STACKSIZE 50000m

"""

if len(sys.argv) > 2:
	np = sys.argv[2]
else:
	np = '1'

if len(sys.argv) > 3:
        mod = sys.argv[3]
else:
        mod = 'medium-mem-1-s'

st = st.replace('file', file)
st = st.replace('proc', np)
st = st.replace('mod', mod)

st1 = '/users/nirajv/xtb_exe/bin/xtb ' +file + ' --opt tight\n'

f.write(st+st1) 
f.close()
