#! /usr/bin/env python

import os
import sys
import screed
from collections import defaultdict

log = str(snakemake.log)
fasta = str(snakemake.input.fasta)
og2genes = str(snakemake.input.og2genes)

default_prefix = os.path.basename(fasta).rsplit('.fa',1)[0]

prefix = snakemake.params.get("prefix", default_prefix)

##### fix these
output_dirname = str(snakemake.output.fastadir)
sketch_dirname = str(snakemake.output.sketchdir)

if not os.path.exists(output_dirname):
    try:
        os.makedirs(output_dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


# some seqs in mult groups? use dict
gene2ogD = defaultdict(list)
with open(og2genes, "r") as ogInfo:
    for line in ogInfo:
        og, gene = line.strip().split("\t")
        gene2ogD[gene].append(og)


with open(log, "w") as out_log:
    out_log.write(f"Splitting {fasta} by orthoDB gene group. Writing files to {output_dirname} \n")

    ogD = defaultdict(list)

    for n, record in enumerate(screed.open(fasta)):
        if n > 0 and n % 100000 == 0:
            out_log.write(f"working on {str(n)}th contig\n")
        #this_og = record.name.rsplit("\t", 1)[1]
        this_gene = record.name.rsplit("\t", 1)[0]
        # get orthogroups this is a part of
        ogs = gene2ogD[this_gene]
        for og in ogs:
           # add gene to each orthogroup
           ogD[og].append(record)

        #this_gene = og2geneD[this_og]
        #ogD[this_og].append(record)
        #ogD[this_gene].append(record)

    filenum = 0
    for orthogroup, og_records in ogD.items():
        outfile = os.path.join(output_dirname, f"{prefix}.{orthogroup}.fa")
        filenum+=1
        with open(outfile, "w") as out:
            for record in og_records:
                out.write(f">{record.name}\n{record.sequence}\n")

        # now sketch this file
        #sketchfile = os.path.join(sketch_dirname, f"{prefix}.{orthogroup}.sig")


    out_log.write(f"{str(n)} contigs written to {str(filenum)} individual orthogroup fasta files\n")
