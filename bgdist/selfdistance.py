import statistics
import sys

from grimoire.genome import Reader

genome = Reader(gff=sys.argv[2], fasta=sys.argv[1])
ratios = []
disttot = 0
gcount = 0
m2_dtot = 0
for chrom in genome:
	#print(chrom.name, file=sys.stderr, flush=True)
	icount = {}
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		icount[(f.beg, f.end)] = f.score

	for gene in chrom.ftable.build_genes():
		if len(gene.transcripts()) == 0: continue # skip ncRNAs
		if gene.issues: continue # skip genes with obvious oddities

		for tx in gene.transcripts():
			bad_tx = False
			counts = []
			for intron in tx.introns:
				isig = (intron.beg, intron.end)
				if isig not in icount:
					bad_tx = True
					break
				counts.append(icount[isig])

		if len(counts) < 2 or bad_tx:
			continue
		
		if sys.argv[3] == 'TRUE':	
			print(gene.id, counts)
			continue

		# Method 1
		countmean = statistics.mean(counts)

		cratios = [count/countmean  for count in counts]
		distsum = sum([abs(cratio - 1) for cratio in cratios])
		gcount += 1
		disttot += distsum

		# Method 2
		countot = sum(counts)

		cfreqs = [count/countot for count in counts]

		unidist = [1/len(counts) for i in range(len(counts))]
		
		m2_dsum = sum([abs(cfreq - unifreq) for cfreq, unifreq in zip(cfreqs, unidist)])
		m2_dtot += m2_dsum

		

print('Background Distance, Method 1:', disttot/gcount)
print('Background Distance, Method 2:', m2_dtot/gcount)

"""
			if len(counts) < 2 or bad_tx: continue
			for i in range(0, len(counts)):
				for j in range(i+1, len(counts)):
					rmax = max(counts[i]/counts[j], counts[j]/counts[i])
					if rmax > 2: continue # 2-fold max
					ratios.append(rmax)

print(statistics.mean(ratios), statistics.stdev(ratios), statistics.median(ratios))
#for ratio in ratios:
#	print(f'{ratio:.2f}')
"""