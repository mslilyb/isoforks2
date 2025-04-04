import argparse
import os



"""
Script that takes .gff3 files from a target directory and outputs field of choice.


"""

parser = argparse.ArgumentParser(description=f'Script that pulls out terms of interest from WormBase .gff3 files in a directory.',
	epilog= 'The default behavior is to search the current directory for WormBase gene IDs. Other terms can be specified')


# Options

parser.add_argument('-d','--directory', required=False, type=str, metavar='<file>',default=None,
	help='Directory containing .gff3 files to pull from. Default: Current Directory')
parser.add_argument('-t','--target', required=False, type=str, metavar='<field>', default='gene',
	help='Name of target field desired to match. Default: %(default)s')

parser.add_argument('--ID', required=False, type=str, metavar='<field>', default=None,
	help='ID of desired information from matching line. Defaults to printing entire line.')



args = parser.parse_args()
path = None
if not args.directory:
	path = os.getcwd()
else:
	path = os.path.abspath(args.directory)


gffs = [file.path for file in os.scandir(path) if file.path.endswith('.gff3')]

for gff in gffs:
	with open(gff) as fp:
		for line in fp:
			fields = line.rstrip().split()
			out = line
			if args.target in line:
				if args.ID:
					for field in fields:
						if args.ID in field:
							out = field

				print(out)

