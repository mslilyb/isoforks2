import argparse
import sys

import isoform2
from grimoire.genome import Reader

parser = argparse.ArgumentParser()
parser.add_argument('model', help='splicing model file (from smallgenes)')
parser.add_argument('fasta', help='fasta file (from smallgenes)')
parser.add_argument('gff', help='gff file (from smallgenes)')
parser.add_argument('--min-orf', type=int, default=25,
	help='minimum distance from start to stop [%(default)i]')
parser.add_argument('--nostop', type=float, default=1e-3,
	help='degredation from no stop codon [%(default)g]')
parser.add_argument('--ejc', type=float, default=1e-3,
	help='degredation from ejc after stop codon [%(default)g]')
parser.add_argument('--tail', type=float, default=1e-3)
#parser.add_argument('--kozak', type=int, default=
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

print(f'CDS: {cds_beg}..{cds_end}') # 114..143, 190..351

for isoform in locus.isoforms:
	#print(isoform.score)
	cds_found = False
	cds_idx = None
	for i, (beg, end) in enumerate(isoform.exons):
		if cds_beg >= beg and cds_beg <= end:
			cds_found = True
			cds_idx = i			
	if not cds_found:  # use alt start based on longest ORF?
		continue       # for now, the transcript won't translate

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
		#print('i', exon_beg, exon_end)
		cds_seqs.append(region.seq[exon_beg:exon_end+1])
	cds_seq = ''.join(cds_seqs)
	protein = isoform2.translate_str(cds_seq)
	stop = protein.find('*')
	print(protein, cds_seqs)
	if stop == -1: # this is a target for non-stop degredation
		print('no stop codon')
		continue   # change this at some point
	orf_len = stop * 3 + 2
	stop_pos = cds_beg + orf_len + intron_len
	if stop_pos < cds_end: print('NMD target')
	elif stop_pos > cds_end: print('normal stop skipped')
	else: print('expected')
	

	# is the nomral cds_beg inside the transcript?
	# if so, use it. if not???


#name, seq = next(isoform2.read_fasta(arg.fasta))
#gene = isoform2.Locus(name, seq, model)
# transcripts without stop codons would also be destroyed



sys.exit('dev')

isoforms = []
for isoform in gene.isoforms:
	print(isoform.beg)
	