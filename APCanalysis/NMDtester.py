import argparse
import re
import sys

import isoform2
from grimoire.genome import Reader

##### Copy-paste from Ian's nmd-ish.py #####
# need to get gff file from Locus
# re-order gff with new probabilities
# ch.2_1 has lots of nmd targets and non stops for testing at limit 10

parser = argparse.ArgumentParser(description='Does nmd-ish improve APC predictions?')
parser.add_argument('model', help='splice model file')
parser.add_argument('fasta', help='fasta file')
parser.add_argument('gff', help='gff file')
parser.add_argument('--min-orf', type=int, default=25,
	help='minimum distance from start to stop [%(default)i]')
parser.add_argument('--nonstop', type=float, default=1e-3,
	help='degredation from no stop codon [%(default)g]')
parser.add_argument('--nmd', type=float, default=1e-3,
	help='degredation from NMD [%(default)g]')
parser.add_argument('--min-ejc', type=int, default=10,
	help='NMD minimum distance from stop to ejc [%(default)i]')
parser.add_argument('--utr', type=int, default=300,
	help='NMD long tail trigger length [%(default)i]')
parser.add_argument('--limit', type=int, default=100,
	help='maximum number of isoforms [%(default)i]')
parser.add_argument('--flank', type=int, default=99,
	help='flanking, non-coding sequence [%(default)i]')

args = parser.parse_args()

model = isoform2.read_splicemodel(args.model)
reader = Reader(fasta=args.fasta, gff=args.gff)
region = next(reader)
locus = isoform2.Locus(region.name, region.seq, model, limit=args.limit)

# create gff to rewrite it
# don't want to edit code in class Locus...
with open('locus.tmp', 'w') as tmp:
	locus.write_gff(tmp)

headers = {}
isos = {}
count = 0
with open('locus.tmp', 'r') as gfp:
	for line in gfp.readlines():
		line = line.rstrip()
		if line.startswith('#'): 
			hline = line.split(' ')
			headers[hline[1][:-1]] = hline[2]
			continue
		if line == '': 
			count += 1
			continue
		if count not in isos:
			isos[count] = []
			isos[count].append(line)
		else:
			isos[count].append(line)
			
print(headers)

for i in isos:
	print(i)
	for j in isos[i]:
		print(j)




'''
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
atgs = set()
for tx in txs:
	cdss = sorted(tx.cdss, key=lambda x: x.beg)
	atgs.add(cdss[0].beg -1)
cds_beg = sorted(list(atgs))[0]

prevs = [iso.prob for iso in locus.isoforms]
posts = []
rtypes = []
for iso in locus.isoforms:
	iso.translate(cds_beg)
	if   iso.rnatype == 'non-stop': iso.prob *= 1e-3
	elif iso.rnatype == 'nmd-target': iso.prob *= 1e-3
	posts.append(iso.prob)
	rtypes.append(iso.rnatype)
	
for a, b in zip(prevs, posts):
	print(a, b)
	
# prevs is all the probabilites before labeling
# posts is all the probabilites after labeling
# labeling with non-stop or nmd-target multiple the score by 1e-3

# the probability in iso.prob has been changed if non-stop or nmd-target
#for iso in locus.isoforms:
#	print(iso.prob)

total = sum([iso.prob for iso in locus.isoforms])
for iso in locus.isoforms:
	iso.prob /= total
	
for iso in locus.isoforms:
	print(iso.prob)
print('##########')
print('rna-type', 'original', 'reduced', 'normalized', sep='\t')
for rtype, prev, post, iso in zip(rtypes, prevs, posts, locus.isoforms):
	print(f'{rtype}\t{prev:.3g}\t{post:.3g}\t{iso.prob:.3g}')

# now i need to recompute the Mdist using the new probabilities
# start without weighted models
# need to read in nmd-ish output
'''

'''
python3 nmd-ish.py models/worm.splicemodel ../datacore2024/project_splicing/smallgenes/ch.2_1.fa ../datacore2024/project_splicing/smallgenes/ch.2_1.gff3 --limit 10
'''

parser = argparse.ArgumentParser()
parser.add_argument('nmd_file', help='nmd-ish output file')

