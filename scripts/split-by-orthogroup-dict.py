#! /usr/bin/env python

import os
import sys
import screed
from collections import defaultdict

log = str(snakemake.log)
fasta = str(snakemake.input)

default_prefix = os.path.basename(fasta).rsplit('.fa',1)[0]

prefix = snakemake.params.get("prefix", default_prefix)
output_dirname = str(snakemake.output)

if not os.path.exists(output_dirname):
    try:
        os.makedirs(output_dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

with open(log, "w") as out_log:
    out_log.write(f"Splitting {fasta} by orthogroup. Writing files to {output_dirname} \n")

    ogD = defaultdict(list)

    for n, record in enumerate(screed.open(fasta)):
        if n > 0 and n % 100000 == 0:
            out_log.write(f"working on {str(n)}th contig\n")
        this_og = record.name.rsplit("\t", 1)[1]
        ogD[this_og].append(record)

    filenum = 0
    for orthogroup, og_records in ogD.items():
        outfile = os.path.join(output_dirname, f"{prefix}-{orthogroup}.fa")
        filenum+=1
        with open(outfile, "w") as out:
            for record in og_records:
                out.write(f">{record.name}\n{record.sequence}\n")

    out_log.write(f"{str(n)} contigs written to {str(filenum)} individual orthogroup fasta files\n")
