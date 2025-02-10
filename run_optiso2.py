import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='wrapper for optiso2')
parser.add_argument('apc_dir', type=str, metavar='<directory>', 
    help='directory with APC .gff3 and .fa files')
parser.add_argument('model', type=str, metavar='<file>',
    help='splice model file')
parser.add_argument('--program', required=False, type=str, default='optiso2', 
    metavar='<executable>', help='path to optiso2')

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

cmds = []
for gID in fpaths:
    cmd = f'./{args.program} {fpaths[gID][0]} {fpaths[gID][1]} {args.model}'
    cmds.append(cmd)

with open('results_optiso2_defaults.txt', 'w') as file:
    for c in cmds:
        c = c.split(' ')
        gid = c[1].split('.')[-2]
        print(f'working on {gid}...')
        result = subprocess.run(c, stdout=subprocess.PIPE, text=True)
        print(result)


