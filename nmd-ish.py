import argparse
import re
import sys

import isoform2
from grimoire.genome import Reader

def display_isoform(gseq, tx, cds_beg, wrap=100, flank=99):
	
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
parser.add_argument('--flank', type=int, default=99,
	help='flanking, non-coding sequence [%(default)i]')
arg = parser.parse_args()

# create the locus
model = isoform2.read_splicemodel(arg.model)
reader = Reader(fasta=arg.fasta, gff=arg.gff)
region = next(reader)
locus = isoform2.Locus(region.name, region.seq, model)

# find the canonical start codon (5' atg if there are more than one)
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
atgs = set()
for tx in txs:
	cdss = sorted(tx.cdss, key=lambda x: x.beg)
	atgs.add(cdss[0].beg -1)
cds_beg = sorted(list(atgs))[0]

# examine the isoforms to determine if they are NMD targets
for isoform in locus.isoforms:
	
	# prefer the annotated atg to the first atg
	atg_found = True
	atg_pos = None
	for i in range(cds_beg, cds_beg + 3):
		if isoform.dnaidx_to_rnaidx(i) is None:
			atg_found = False
			break
	if atg_found is False: # find the first atg of the transcript
		x = isoform.txseq.find('ATG')
		if x != -1: atg_pos = isoform.rnaidx_to_dnaidx(x)
	else: atg_pos = cds_beg
	
	# short-circuit on non-start transcripts (maybe not be degraded)
	if atg_pos is None: continue
	
	# check for NMD and other surveillance
	isoform.translate(atg_pos)

	print(isoform.start, isoform.stop)

	display_isoform(region.seq, isoform, atg_pos)
	
	continue


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