import argparse
import isoform2
import os
import subprocess
import sys

apcgffs = os.listdir(os.path.abspath(sys.argv[1]))

rnagffs = []
for file in os.listdir(os.path.abspath(sys.argv[2])):
	if file.endswith('.gff3'):
		rnagffs.append(file)



apcgs = [os.path.abspath(sys.argv[1]) + '/' + file for file in apcgffs]
rnags = [os.path.abspath(sys.argv[2]) + '/' + file for file in rnagffs]
apcgs.sort()
rnags.sort()
apcgffs.sort()

assert(len(apcgs) == len(rnags))


count = 0
misscount = 0

for f1, f2, f3 in zip(apcgs, rnags, apcgffs):
	misscount = 0
	apcis = isoform2.get_introns(os.path.abspath(f1))
	rnais = isoform2.get_introns(os.path.abspath(f2))

	name = f3.rstrip('.APC.gff3')
	apcfreqs = []
	rnafreqs = []

	for key in rnais.keys():
		rnafreqs.append(rnais[key])
		if key not in apcis.keys():
			misscount += 1
			apcfreqs.append(0.0)

		else:
			apcfreqs.append(apcis[key])

	assert(len(rnafreqs) == len(apcfreqs))

	ustat, pval = isoform2.mannu(apcfreqs, rnafreqs)

	print(name, ustat, pval, sep=',')
	print('missed', misscount/len(rnais.keys()) * 100, 'percent of introns', file=sys.stderr)




'''
Needs documentation. Here was the old method for comparison:

def mannu(p, q):
	"""wrapper for mannuwhitney, asserts included. unknown if needed"""
	newp = [i for i in p if not math.isclose(i, 0, abs_tol=1e-6)]
	newq = [i for i in q if not math.isclose(i, 0, abs_tol=1e-6)]
	return scistats.mannwhitneyu(newp, newq, alternative='two-sided')

this involved dropping 0s from BOTH sets of freqs, which explains the difference in results. i believe we should add 0s to the apc data to acknowledge missed introns
'''



