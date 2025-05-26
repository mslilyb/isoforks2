import argparse
import isoform2
import os
from scipy import stats
import subprocess
import sys

parser = argparse.ArgumentParser(description=f'Tool to test various statistics using Kolgorov-Smirnov tests')
parser.add_argument('apcd', type=str, metavar='<directory>', help='path to APC .gffs')
parser.add_argument('sgd', type=str, metavar='<directory>', help='path to smallgenes.gffs')

# Flags
parser.add_argument('-a', '--adjust', action='store_true', help='account for normality in frequencies \
	force enables the -n flag.')
parser.add_argument('-n','--normal', action='store_true', help='test to see if either set of \
	intron frequencies could reasonbly come from a normal distribution')

args = parser.parse_args()

apcgffs = os.listdir(os.path.abspath(args.apcd))

sggffs = []
for file in os.listdir(os.path.abspath(args.sgd)):
	if file.endswith('.gff3'):
		sggffs.append(file)



apcgs = [os.path.abspath(args.apcd) + '/' + file for file in apcgffs]
sggs = [os.path.abspath(args.sgd) + '/' + file for file in sggffs]
apcgs.sort()
sggs.sort()
apcgffs.sort()

assert(len(apcgs) == len(sggs))


count = 0
misscount = 0
smallcount = 0
skipped = []

for f1, f2, f3 in zip(apcgs, sggs, apcgffs):
	missedit = False
	apcis = isoform2.get_introns(os.path.abspath(f1))
	sgis = isoform2.get_introns(os.path.abspath(f2))

	apcfreqs = []
	sgfreqs = []
	apc_isnorm = False
	sg_isnorm = False

	for key in sgis.keys():
		sgfreqs.append(sgis[key])
		if key not in apcis.keys():
			missedit = True
			apcis[key] = 0.0
		
		apcfreqs.append(apcis[key])

	assert(len(sgfreqs) == len(apcfreqs))

	if missedit:
		misscount += 1


	if args.adjust or args.normal:
		apc_isnorm = stats.kstest(apcfreqs, stats.norm.cdf).pvalue >= 0.05
		sg_isnorm = stats.kstest(sgfreqs, stats.norm.cdf).pvalue >= 0.05

	if apc_isnorm and sg_isnorm and args.adjust:
		stat, pval = stats.ttest_ind(apcfreqs, sgfreqs, nan_policy='raise')
	else:
		stat, pval = stats.kstest(apcfreqs, sgfreqs)

	print(f3, stat, pval, sep=',', end='')

	if args.adjust:
		if apc_isnorm and sg_isnorm:
			print('', 'T-TEST', sep=',')

		else:
			print('','KS-TEST', sep=',')





#print(misscount, "genes where APC failed to predict introns found in WB", file=sys.stderr)

'''
Needs documentation. Here was the old method for comparison:

def mannu(p, q):
	"""wrapper for mannuwhitney, asserts included. unknown if needed"""
	newp = [i for i in p if not math.isclose(i, 0, abs_tol=1e-6)]
	newq = [i for i in q if not math.isclose(i, 0, abs_tol=1e-6)]
	return scistats.mannwhitneyu(newp, newq, altesgtive='two-sided')

this involved dropping 0s from BOTH sets of freqs, which explains the difference in results. i believe we should add 0s to the apc data to acknowledge missed introns
'''



