# code that reads gff file contents into some dictionaries

headers = {}
isos = {}
count = 0
with open('prenmdish.gff', 'r') as gfp:
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