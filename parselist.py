import os

genes = {} #key is gene, value is dict of booleans

for file in os.listdir(os.getcwd()):
	if file.endswith('_genelist'):
		with open(file) as fp:
			used = 'True'
			type = file.split('_')[0]
			for line in fp:
				if 'Included:' in line:
					continue
				elif 'Excluded:' in line:
					used = 'False'
					continue

				gene = line.split('-')[0].rstrip()
				
				if gene not in genes.keys():
					genes[gene] = {'TEA': 'N/A', 'PEA': 'N/A', 'GO': 'N/A'}

				genes[gene][type] = used
				

print(f'genes,TEA,PEA,GO')

for key, values in genes.items():
	outstr = key + ',' + ','.join(values.values())
	print(outstr)

