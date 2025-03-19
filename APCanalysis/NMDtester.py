import argparse
import glob

import isoform2
from grimoire.genome import Reader

##### Copy-paste from Ian's nmd-ish.py #####
# need to get gff file from Locus
# order of isos in gffs does not matter?
# ch.2_1 has lots of nmd targets and non stops for testing at limit 10

parser = argparse.ArgumentParser(description='Does nmd-ish improve APC predictions?')
parser.add_argument('model', help='splice model file')
#parser.add_argument('fasta', help='fasta file')
#parser.add_argument('gff', help='gff file')
parser.add_argument('smallgenes', help='smallgenes directory')
parser.add_argument('--weights', required=False, type=str, 
	help='optiso2 output with model weights')
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



# can I implement the new distance equation?
# Dtf(P||Q) = Sum |pi - qi| * max(pi/qi, qi/pi)
# how to account for 0s?

def compare_dists(model, infasta, ingff, inlimit, inweights=None):

	model = isoform2.read_splicemodel(model)
	reader = Reader(fasta=infasta, gff=ingff)
	region = next(reader)
	if inweights:
		locus = isoform2.Locus(region.name, region.seq, model, limit=inlimit, 
						weights=inweights)
	else:
		locus = isoform2.Locus(region.name, region.seq, model, limit=inlimit)

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
	prevs = [iso.prob for iso in locus.isoforms] # don't need to print here
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
	i1 = isoform2.get_introns(ingff)
	i2 = isoform2.get_introns('prenmdish.gff.tmp')
	i3 = isoform2.get_introns('postnmdish.gff.tmp')

	dist1, details1 = isoform2.expdiff(i1, i2)
	dist2, details2 = isoform2.expdiff(i1, i3)

	info = [region.name, dist1, dist2, dist2-dist1]

	return info

file_pairs = {}
for fpath in glob.glob(f'{args.smallgenes}/*.gff3'):
	gid = fpath.split('/')[-1].split('.')[1]
	file_pairs[gid] = [fpath]

for fpath in glob.glob(f'{args.smallgenes}/*.fa'):
	gid = fpath.split('/')[-1].split('.')[1]
	file_pairs[gid].append(fpath)

if args.weights:
	weights = {}
	with open(args.weights, 'r') as file:
		for line in file.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split(',')
			weights[line[0]] = [x for x in line[1:]]
	
	print(f'gene\tprenmdish\tpostnmdish\tdelta')
	for p in file_pairs:
		w = weights[p]
		wf = {
			'wacc': float(w[1]),
			'wdon': float(w[2]),
			'wexs': float(w[3]),
			'wins': float(w[4]),
			'wexl': float(w[5]),
			'winl': float(w[6]),
			'winf': float(w[7])
		}
		info = compare_dists(args.model, file_pairs[p][1], file_pairs[p][0],
					args.limit, inweights=wf)
		print(f'{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}')

else:
	print(f'gene\tprenmdish\tpostnmdish\tdelta')
	for p in file_pairs:
		info = compare_dists(args.model, file_pairs[p][1], file_pairs[p][0],
					args.limit)
		print(f'{info[0]}\t{info[1]}\t{info[2]}\t{info[3]}')


# why are the Mdists worse after including the weights???
# need to parellize


