import argparse

parser = argparse.ArgumentParser()
parser.add_argument('genome')
parser.add_argument('gff')

args = parser.parse_args()

introns = []
exons = []
with open(args.gff, 'r') as gp:
	for line in gp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if line[2] == 'intron':
			intron = (line[0], int(line[3]), int(line[4]), line[6])
			#if line[6] == '-' and int(line[4])-int(line[3]) < 100:
				#print(line)
			introns.append(intron)
		if line[2] == 'exon':
			exon = (line[0], int(line[3]), int(line[4]))
			exons.append(exon)

chroms = {}
with open(args.genome, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			chrom = line.split(' ')[0][1:]
		else:
			if chrom not in chroms:
				chroms[chrom] = line
			else:
				chroms[chrom] += line

accs = []
dons = []
eseqs = []
iseqs = []
exonsum = 0

total = 0
nc = 0
for i in introns:
	beg = i[1]-1
	end = i[2]
	iseq = chroms[i[0]][beg:end]
	aseq = iseq[-6:]
	dseq = iseq[:5]
	if aseq.endswith('AG'): total += 1
	if not aseq.endswith('AG'):
		#if i[2] - i[1] < 100:
            #print(i, aseq)
		nc += 1
	accs.append(aseq)
	dons.append(dseq)
	iseqs.append((iseq, i[3]))

print(total)
print(nc)

compseq = ''

for i in iseqs:
	if i[1] == '-':
		compseq = i[0]
		break
		
print(compseq)
	
	
# need to weight the introns, so highly expressed ones are represented more
'''
for e in exons:
	beg = e[1]-1
	end = e[2]
	eseq = chroms[i[0]][beg:end]
	exonsum += len(eseq)
	eseqs.append(eseq)
'''
