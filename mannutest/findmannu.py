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
smallcount = 0
skipped = []

for f1, f2, f3 in zip(apcgs, rnags, apcgffs):
	misscount = 0
	apcis = isoform2.get_introns(os.path.abspath(f1))
	rnais = isoform2.get_introns(os.path.abspath(f2))

	apcfreqs = []
	rnafreqs = []

	for key in rnais.keys():
		rnafreqs.append(rnais[key])
		if key not in apcis.keys():
			misscount += 1
			apcis[key] = 0.0
			apcfreqs.append(apcis[key])
		else:
			apcfreqs.append(apcis[key])

	assert(len(rnafreqs) == len(apcfreqs))

	if len(rnafreqs) * len(apcfreqs) < 16:
		smallcount += 1
		for key in rnais.keys():
			print(f1, len(rnafreqs), key, rnais[key], apcis[key], rnais[key] - apcis[key], sep=',', file=sys.stderr)
		continue #failed to meet minimum sample size for mann U 2 tailed @ 5%

	ustat, pval = isoform2.mannu(apcfreqs, rnafreqs)

	print(f1, ustat, pval, sep=',')


print(smallcount, "skipped genes with insufficient sample sizes.", file=sys.stderr)

'''
Needs documentation. Here was the old method for comparison:

def mannu(p, q):
	"""wrapper for mannuwhitney, asserts included. unknown if needed"""
	newp = [i for i in p if not math.isclose(i, 0, abs_tol=1e-6)]
	newq = [i for i in q if not math.isclose(i, 0, abs_tol=1e-6)]
	return scistats.mannwhitneyu(newp, newq, alternative='two-sided')

this involved dropping 0s from BOTH sets of freqs, which explains the difference in results. i believe we should add 0s to the apc data to acknowledge missed introns
'''



