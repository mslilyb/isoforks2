import sys

introns = sys.argv[1]

minin = 35

pds = []
pas = []
nt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
total = 0
with open(introns, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		pd = line[:17]
		if pd.startswith('GT'):
			pds.append(pd)
		pa = line[-18:]
		if pa.endswith('AG'):
			pas.append(pa)
		for nt in line:
			nt_counts[nt] += 1
			total += 1	
	
f = [nt_counts[x]/total for x in nt_counts]
nseqs = len(pds)
bg = [nseqs * f[0], nseqs * f[1], nseqs * f[2], nseqs * f[3]]
print(f)

# instead just make a table that can be used in R	
'''
def chi_test(seqs, background):
	
	count = []
	total = 0
	for seq in seqs:
		print(seq, '#')
		
chi_test(pds, bg)
'''
		
	
count = []
total = 0
for seq in pds:
	for i, nt in enumerate(seq):
		if len(count) <= i:
			count.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
		count[i][nt] += 1
		total += 1

print(count)

chis = []
for site in count:
	chi = 0
	for i, nt in enumerate(site):
		# even if GT is 100%, background for G and T is different
		c = ((site[nt] - bg[i]) ** 2) / bg[i]
		chi += c
	chis.append(chi)
	
for x in chis:
	print(x)
		


'''
seqs = [
	'GTAAG',
	'GTACG',
	'GTCCG',
	'GTAAG',
	'GTCAC'
]

count = []
for seq in seqs:
	for i, nt in enumerate(seq):
		if len(count) <= i:
			count.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
		count[i][nt] += 1
	
print(count)

# chi square test
counts = {'A': 83, 'C': 30, 'G': 49, 'T': 84}
freqs = {}
total = 0
for n in counts:
	total += counts[n]
	
print(total)

for n in counts:
	if n not in freqs:
		freqs[n] = counts[n]/total
		
print(freqs)
		
bg = {'A': 0.3279, 'C': 0.1915, 'G': 0.2043, 'T': 0.2763}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
print(ecounts)
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# chi square 2
counts = {'A': 103, 'C': 44, 'G': 46, 'T': 53}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# chi square 2
counts = {'A': 121, 'C': 36, 'G': 38, 'T': 51}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# alpha for rejection
a = 0.003012705

# site 1 P is 0.0177
# 3 degrees of freedom? 4 - 1 = 3
'''



















