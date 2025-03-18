import argparse
import re
import sys

import isoform2
from grimoire.genome import Reader

def display_isoform(gseq, tx, cds_beg, wrap=80, flank=99):
	# used for debugging
	print(iso.aaseq)
	rna = [' '] * len(gseq)
	for beg, end in tx.exons:
		for i in range(beg+1, end): rna[i] = 'X'
		rna[beg] = '['
		rna[end] = ']'
	for beg, end in tx.introns:
		for i in range(beg, end+1): rna[i] = '-'
	tseq = ''.join(rna)
	txcoor = 0
	tpos = []
	for nt in tseq:
		if nt != '-': txcoor += 1
		tpos.append(txcoor)

	cds = [' '] * len(gseq)
	x = cds_beg - tx.exons[0][0] # first atg
	while True:
		c0 = tx.rnaidx_to_dnaidx(x)
		c1 = tx.rnaidx_to_dnaidx(x+1)
		c2 = tx.rnaidx_to_dnaidx(x+2)
		if c0 is None: break
		if c1 is None: break
		if c2 is None: break
		codon = gseq[c0] + gseq[c1] + gseq[c2]
		aa = isoform2.GCODE[codon]
		cds[c1] = aa
		if codon == 'TAA' or codon == 'TAG' or codon == 'TGA': break
		x += 3
	cseq = ''.join(cds)

	for i in range(0, len(gseq), wrap):
		rbeg = i
		rend = i + wrap

		# genome coor
		for j in range(0, wrap, 10):
			print(f'{i+j+10:>10}', end='')
		print()
		for j in range(0, wrap, 10):
			print(' ' * 9, end='')
			print('|', end='')
		print()

		# sequences
		print(gseq[i:i+wrap])
		has_cds = False
		has_exon = False
		if re.search(r'\S', cseq[i:i+wrap]): has_cds = True
		if re.search(r'\S', tseq[i:i+wrap]): has_exon = True

		if has_cds: print(cseq[i:i+wrap])

		if has_exon:
			print(tseq[i:i+wrap])
			for j in range(0, wrap, 10):
				print(' ' * 9, end='')
				print('|', end='')
			print()

			for j in range(0, wrap, 10):
				p = i + j + 9
				if p >= len(tseq): break
				x = tpos[p]

				if tseq[p] == '.': print(' ' * 10, end='')
				else: print(f'{x:>10}', end='')
			print()
		print()
	sys.exit()

parser = argparse.ArgumentParser(description='Does nmd-ish improve APC predictions?')
parser.add_argument('model', help='splice model file')
parser.add_argument('fasta', help='fasta file')
parser.add_argument('gff', help='gff file')
parser.add_argument('--debug', action='store_true')

args = parser.parse_args()

# ch.2_1 has lots of nmd targets and non stops

# what does nmd-ish do? let's recreate it

model = isoform2.read_splicemodel(args.model)
reader = Reader(fasta=args.fasta, gff=args.gff)
region = next(reader)
locus = isoform2.Locus(region.name, region.seq, model, limit=10)

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
	print(iso.rnatype, iso.prob)
	if args.debug: display_isoform(region.seq, iso, cds_beg)
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