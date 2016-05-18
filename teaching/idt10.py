#!/usr/bin/env python

from Bio import Entrez, SwissProt
Entrez.email="name@domain.com"

# get info about NCBI databases
handle = Entrez.einfo()
record = Entrez.read(handle)
print record

# define database
db="protein"
# look for `opsin 1` from human
query = '"tp63" AND "homo sapiens"[organism]'
handle = Entrez.esearch(db=db, retmax=10, term=query)
record = Entrez.read(handle)
print record

# get genbank ID associated with above query
gid = record['IdList'][0]

# fetch protein sequence in FASTA
handle = Entrez.efetch(db=db, id=gid, rettype="fasta")
fasta = "".join(handle.readlines())
print fasta


###
# perform blastP search remotely
from Bio.Blast import NCBIWWW, NCBIXML

myEntrezQuery = "Homo sapiens[Organism] OR Mus musculus[Organism] OR Gallus gallus[Organism] OR Danio rerio[Organism] OR Drosophila melanogaster[Organism]"
result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta, expect=1e-25, hitlist_size=100, entrez_query=myEntrezQuery)
blast_record = NCBIXML.read(result_handle)
# dt=datetime.now()-t0

# you can explore blast_record object using `blast_record.` + TAB

# parse alignemnts
for alignment in blast_record.alignments:
  for hsp in alignment.hsps:
    print alignment.title, hsp.expect
  
blast_hsp = blast_qresult[0][0]
print blast_hsp
print(blast_hsp.aln)

# fetch sequences
from Bio import SeqIO
# get gi ids
gids = [a.hit_id.split('|')[1] for a in blast_record.alignments]

handle = Entrez.efetch(db=db, id=",".join(gids), rettype="fasta")
records = []
for r in SeqIO.parse(handle, "fasta"):
  r.id = r.id.split('|')[-1]
  #r.description = ""
  records.append(r)

# MSA using muscle
import subprocess, sys
from Bio.Align.Applications import MuscleCommandline

fn="seqs.fa"
with open(fn, "w") as out:
  SeqIO.write(records, out, "fasta")
  
muscle_cmd = MuscleCommandline(clwstrict=True, input=fn)
print(muscle_cmd)
child = subprocess.Popen(str(muscle_cmd), stdout=subprocess.PIPE, shell=(sys.platform!="win32"))
child.wait()

# read alignments
from Bio import AlignIO
align = AlignIO.read(child.stdout, "clustal")
print(align)

# convert into PHYLIP format
phylip="seqs.phy"
with open(phylip, 'w') as out:
  AlignIO.write(align, out, 'phylip')

###
# reconstruct phylogenetic tree
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline, FastTreeCommandline
#cmd = PhymlCommandline(input=phylip, datatype='aa')
cmd = FastTreeCommandline(input=phylip, out=phylip+".nw")
out_log, err_log = cmd()

tree = Phylo.read(phylip+".nw", 'newick') # '_phyml_tree.txt'
Phylo.draw_ascii(tree)

# ete
import ete3
t=ete3.PhyloTree(phylip+".nw")

# root by mid-point
t.set_outgroup(t.get_midpoint_outgroup())

print(t)
t.show()
t.render(fn+".svg")


