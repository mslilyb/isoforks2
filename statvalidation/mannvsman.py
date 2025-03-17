import sys
import os

resultfiles = os.listdir(os.path.abspath(sys.argv[1]))
print(f'gene,manhattan,weird')
for rfile in resultfiles:
	fname = rfile.rstrip('.cmp')
	with open(f'{sys.argv[1]}/{rfile}') as fp:
		resul = []
		almost_done = False
		for line in fp:

			if line.startswith('mann-u'):
				almost_done = True

			if line.startswith('test'):
				slots = line.split(f'\t')
				resul.append(slots[-1].rstrip())

			elif line.startswith('p-value'):
				slots = line.split(f'\t')
				resul.append(slots[-1].rstrip())
				if almost_done:
					break
			elif line.startswith('weird:'):
				slots = line.split()
				resul.append(slots[-1].rstrip())
				break
			elif '.' in line:
				manhattan = line.rstrip()

	print(f'{fname},{manhattan},{resul[0]}')



