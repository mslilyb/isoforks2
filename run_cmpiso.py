import os
import subprocess
import sys

apcgffs = os.listdir(os.path.abspath(sys.argv[1]))

rnagffs = []
for file in os.listdir(os.path.abspath(sys.argv[2])):
	if file.endswith('.gff3'):
		rnagffs.append(file)

outdir = sys.argv[3]

os.makedirs(f'{outdir}/', exist_ok=True)
apcgffs.sort()
rnagffs.sort()

test = 5

count = 0

for f1, f2 in zip(apcgffs, rnagffs):
	name = f1.rstrip('.APC.gff3')
	cmd = f'./cmpiso smallgenes/{f2} APCisos/{f1} > {outdir}/{name}.cmp'
	print('running:', cmd)
	result = subprocess.run(cmd, shell=True)

	if result == 1:
		print('uh oh')

	#count += 1

	#if count > test:
		#break