import argparse
import re
import sys

import isoform2
from grimoire.genome import Reader

##### Copy-paste from Ian's nmd-ish.py #####
# need to get gff file from Locus
# order of isos in gffs does not matter?
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

with open('prenmdish.gff.tmp', 'w') as fp:
	locus.write_gff(fp)
			
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
atgs = set()
for tx in txs:
	cdss = sorted(tx.cdss, key=lambda x: x.beg)
	atgs.add(cdss[0].beg -1)
cds_beg = sorted(list(atgs))[0]

# prevs is all the probabilites before labeling
# posts is all the probabilites after labeling
# labeling with non-stop or nmd-target multiplies the probability by 1e-3
prevs = [iso.prob for iso in locus.isoforms]
posts = []
rtypes = []
for iso in locus.isoforms:
	iso.translate(cds_beg)
	if   iso.rnatype == 'non-stop': iso.prob *= 1e-3
	elif iso.rnatype == 'nmd-target': iso.prob *= 1e-3
	posts.append(iso.prob)
	rtypes.append(iso.rnatype)

total = sum([iso.prob for iso in locus.isoforms])
for iso in locus.isoforms:
	iso.prob /= total
	
with open('postnmdish.gff.tmp', 'w') as fp:
	locus.write_gff(fp)

# compare mdist before and after nmd-ish
i1 = isoform2.get_introns(args.gff)
i2 = isoform2.get_introns('prenmdish.gff.tmp')
i3 = isoform2.get_introns('postnmdish.gff.tmp')

dist1, details1 = isoform2.expdiff(i1, i2)
dist2, details2 = isoform2.expdiff(i1, i3)

print(region.name, dist1, dist2, dist2-dist1)

'''
python3 nmd-ish.py models/worm.splicemodel ../datacore2024/project_splicing/smallgenes/ch.2_1.fa ../datacore2024/project_splicing/smallgenes/ch.2_1.gff3 --limit 10
'''

