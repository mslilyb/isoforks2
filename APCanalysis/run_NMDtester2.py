import argparse
import os
import csv
import glob
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time

parser = argparse.ArgumentParser(description='wrapper for NMDtester2')
parser.add_argument('smallgenes', type=str, metavar='<directory>', 
    help='directory with APC .gff3 and .fa files')
parser.add_argument('model', type=str, metavar='<file>',
    help='splice model file')
parser.add_argument('--deqn', required=False, type=str, default='dtc',
	metavar='<string>', help='choose which distance equation to use, '
	'i.e. dtc (Manhattan), dtx, dty, dkl (Kullback-Leibler), '
    'or all [%(default)s]')
parser.add_argument('--limit', required=False, type=int, default=100,
	metavar='<int>', help='limit number of transcripts [%(default)i]')
parser.add_argument('--weights', type=str, metavar='<file>',
    help='file with model weights')
parser.add_argument('--program', required=False, type=str, default='NMDtester2', 
    metavar='<executable>', help='path to NMDtester2')

args = parser.parse_args()

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

def generate(
        prog, fasta, model, gff, wacc, wdon, wexs, wins, wexl, winl, winf, 
        limit, deqn
        ):

    cmd = (
        f'./{prog} {model} {fasta} {gff} --wacc {wacc} --wdon {wdon} '
        f'--wexs {wexs} --wins {wins} --wexl {wexl} --winl {winl} '
        f'--winf {winf} --limit {limit} --deqn {deqn}'
    )
    cmd = cmd.split(' ')
    gid = cmd[2].split('.')[-2]
    output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout
    print(f'{gid} done')
    return output.rstrip()

# order of weights
# fitness,wacc,wdon,wexs,wins,wexl,winl,winf
inputs = []
for gene in file_pairs:
    input = [
        args.program, file_pairs[gene][1], args.model, file_pairs[gene][0],
        weights[gene][1], weights[gene][2], weights[gene][3], 
        weights[gene][4], weights[gene][5], weights[gene][6], 
        weights[gene][7], args.limit, args.deqn
    ]
    inputs.append(input)

def worker(input):
    return generate(
        input[0], input[1], input[2], input[3], input[4], input[5], 
        input[6], input[7], input[8], input[9], input[10], input[11],
        input[12]
     )

with Pool(processes=mp.cpu_count()-1) as pool:
    result = pool.map(worker, inputs)

# pre nmd distance, post nmd distance, delta
# if using --deqn all, from top to bottom:
# dtc, dkl, dtx, dty
with open('nmd.res', 'w') as fp:
	for res in result:
		fp.write(f'{res}\n')



