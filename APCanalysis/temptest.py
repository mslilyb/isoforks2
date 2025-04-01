import tempfile
import isoform2

from grimoire.genome import Reader

model = isoform2.read_splicemodel('../models/worm.splicemodel')
reader = Reader(fasta='testgenes/ch.1_101.fa', gff='testgenes/ch.1_101.gff3')
region = next(reader)
locus = isoform2.Locus(region.name, region.seq, model, limit=100)


with open('test.tmp', 'w') as fp:
    locus.write_gff(fp)

with tempfile.NamedTemporaryFile(delete_on_close=False) as fp:
    locus.write_gff(fp)
