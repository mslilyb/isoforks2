import argparse
import sys

import isoform2
from grimoire.genome import Reader

def display_isoform(gseq, cds_beg, cds_end, exons, wrap=100):
	rna = ['-'] * len(gseq)
	for beg, end in exons:
		for i in range(beg, end+1): rna[i] = gseq[i]
	tseq = ''.join(rna)
	txcoor = 0
	tpos = []
	for nt in tseq:
		if nt != '-': txcoor += 1
		tpos.append(txcoor)
	
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
		print(tseq[i:i+wrap])
		
		# transcriptome coor
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
		
		# print start and stop codons
		if cds_beg >= rbeg and cds_beg <= rend:
			diff = cds_beg - rbeg
			s = ' ' * diff
			print(s, 'ATG', sep='')
		if cds_end >= rbeg and cds_end <= rend:
			diff = cds_end - rbeg -2
			s = ' ' * diff
			stop = tseq[cds_end-2:cds_end+1]
			print(s, stop, sep='')
		print()
		

parser = argparse.ArgumentParser()
parser.add_argument('model', help='splicing model file (from smallgenes)')
parser.add_argument('fasta', help='fasta file (from smallgenes)')
parser.add_argument('gff', help='gff file (from smallgenes)')
parser.add_argument('--min-orf', type=int, default=25,
	help='minimum distance from start to stop [%(default)i]')
parser.add_argument('--nostop', type=float, default=1e-3,
	help='degredation from no stop codon [%(default)g]')
parser.add_argument('--min-ejc', type=int, default=10,
	help='minimum distance from stop to ejc [%(default)i]')
parser.add_argument('--ejc', type=float, default=1e-3,
	help='degredation from ejc after stop codon [%(default)g]')
parser.add_argument('--tail', type=float, default=1e-3,
	help='degredation from a too long 3\'UTR [%(default)g]')
parser.add_argument('--utr', type=int, default=300,
	help='long tail trigger length [%(default)i]')
arg = parser.parse_args()

model = isoform2.read_splicemodel(arg.model)
reader = Reader(fasta=arg.fasta, gff=arg.gff)
region = next(reader)
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
if len(txs) > 1: sys.exit('not sure what to do about that')
tx = txs[0]
cdss = sorted(tx.cdss, key=lambda x: x.beg)
cds_beg = cdss[0].beg -1
cds_end = cdss[-1].end -1

locus = isoform2.Locus(region.name, region.seq, model)

for isoform in locus.isoforms:

	display_isoform(region.seq, cds_beg, cds_end, isoform.exons)

	cds_found = False
	cds_idx = None
	for i, (beg, end) in enumerate(isoform.exons):
		if cds_beg >= beg and cds_beg <= end:
			cds_found = True
			cds_idx = i			
	if not cds_found: continue
		# use some other atg at reduced efficiency?
		# set the isoform score to zero?
		# find the longest ORF?
		# leave isoform score as is?

	cds_seqs = []
	intron_len = 0
	for i in range(len(isoform.exons)):
		exon_beg = isoform.exons[i][0]
		exon_end = isoform.exons[i][1]
		if i >= 1:
			intron_beg = isoform.introns[i-1][0]
			intron_end = isoform.introns[i-1][1]
			intron_len += intron_end - intron_beg + 1
	
		if cds_beg > exon_end: continue
		if cds_beg > exon_beg: exon_beg = cds_beg
		cds_seqs.append(region.seq[exon_beg:exon_end+1])

	
	cds_seq = ''.join(cds_seqs)
	protein = isoform2.translate_str(cds_seq) # debugging
	print(cds_seqs) # debugging
	print(protein)  # debugging
	print(isoform.exons)
	print(isoform.introns)

	# look for first in-frame stop codon
	stop_pos = None
	for i in range(0, len(cds_seq) -2, 3):
		codon = cds_seq[i:i+3]
		if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
			stop_pos= i
			break
	
	# look for ejc downstream of stop codon
	ejc_found = False
	ejcs = []
	prev = 0
	for cds_seq in cds_seqs[1:]:
		ejcs.append(prev + len(cds_seq))
		prev += len(cds_seq)
	for ejc in ejcs:
		if ejc > stop_pos:
			ejc_found = True
			break
	
	print('stop:', stop_pos, ', ecjs:', ejcs)
	
	# measure 3'UTR length
	utr_len = 0 if stop_pos is None else len(cds_seq) - stop_pos 
	
	# finalization
	if ejc_found:
		isoform.score *= arg.ejc
		print('ejc found')
	elif utr_len == 0:
		isoform.score *= arg.nonstop
		print('no stop')
	elif utr_len > arg.utr:
		isoform.score *= arg.tail
		print('long tail')
	else:
		print('ok')



sys.exit('dev')