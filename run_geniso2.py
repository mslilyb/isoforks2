import argparse
import os
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time
import csv

parser = argparse.ArgumentParser(description='wrapper for geniso2')
parser.add_argument('apc_dir', type=str, metavar='<directory>', 
    help='directory with APC .gff3 and .fa files')
parser.add_argument('model', type=str, metavar='<file>',
    help='splice model file')
parser.add_argument('--weights', type=str, metavar='<file>',
    help='file with model weights')
parser.add_argument('--program', required=False, type=str, default='geniso2', 
    metavar='<executable>', help='path to geniso2')

args = parser.parse_args()

fpaths = {}
for file in os.listdir(args.apc_dir):
    gID = file.split('.')[1]
    if gID not in fpaths:
        fpaths[gID] = [(), ()]
    if file.endswith('.fa'):
        fpaths[gID][0] = f'{args.apc_dir}{file}'
    if file.endswith('.gff3'):
        fpaths[gID][1] = f'{args.apc_dir}{file}'

weights = {}
with open(args.weights, 'r') as file:
    for line in file.readlines():
        line = line.rstrip()
        if line.startswith('#'): continue
        line = line.split(',')
        weights[line[0]] = [x for x in line[1:]]


def generate(prog, fasta, model, wdon, wacc, wexs, wins, wexl, winl, winf):

    cmd = (
        f'./{prog} {fasta} {model} --wdon {wdon} --wacc {wacc} '
        f'--wexs {wexs} --wins {wins} --wexl {wexl} --winl {winl} '
        f'--winf {winf}'
    )
    cmd = cmd.split(' ')
    print(cmd)
    gid = cmd[1].split('.')[-2]
    print(gid)
    '''
    output = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout.split()
    print(f'working on {gid}...')
    return [gid] + output
    '''

for f in fpaths:
    #print(f, fpaths[f])
    generate(args.program, fpaths[f][0], args.model, 1, 1, 1, 1, 1, 1, 1)


'''
def worker(inputs):

    return optimize(inputs[0], inputs[1], inputs[2], inputs[3], args.limit)

inputs = []
for gID in fpaths:
    input = [args.program, fpaths[gID][0], fpaths[gID][1], args.model]
    inputs.append(input)

#s = time.perf_counter()
with Pool(processes=mp.cpu_count()-1) as pool:
    result = pool.map(worker, inputs)
    #print(result)
#e = time.perf_counter()
#print('multi:', e-s)
'''