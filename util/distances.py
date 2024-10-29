import random
import math

def freq_dist(vals, pseudo=0):
	assert(pseudo >= 0)
	total = sum(vals) + pseudo * len(vals)
	return [(x+pseudo) / total for x in vals]


def d1(P, Q): # manhattan
	return sum([abs(p - q) for p, q in zip(P, Q)])

def d2(P, Q): # cartesian
	return math.sqrt(sum([(p-q)**2 for p, q in zip(P, Q)]))

def dB(P, Q): # Bhattacharyya
	return -math.log(sum([math.sqrt(p * q) for p, q in zip(P, Q)]))

def dKL(P, Q): # symmetrical version
	return sum([p * math.log2(p/q) + q * math.log2(q/p) for p, q in zip(P, Q)])


data = {
	'd1': [],
	'd2': [],
	'dB': [],
	'dKL': [],
}

for i in range(5000):
	a = []
	b = []
	for _ in range(20):
		a.append(random.randint(0, 20))
		b.append(random.randint(0, 20))
	P = freq_dist(a, pseudo=0.01)
	Q = freq_dist(b, pseudo=0.01)
	data['d1'].append((i, d1(P, Q)))
	data['d2'].append((i, d2(P, Q)))
	data['dB'].append((i, dB(P, Q)))
	data['dKL'].append((i, dKL(P, Q)))
	

for k in data:
	data[k] = sorted(data[k], key=lambda x: x[1])

ks = list(data.keys())
for i in range(len(ks)):
	for j in range(i + 1, len(ks)):
		d = 0
		for x, y in zip(data[ks[i]], data[ks[j]]):
			d += abs(x[0] - y[0])
		print(ks[i], ks[j], d)
		