#! /usr/bin/env python

import os
import sys
import screed

log = str(snakemake.log)
fasta = str(snakemake.input)
contigs_per_file = int(snakemake.params.contigs_per_file)

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
    filenum=0
    out_log.write(f"Splitting {fasta} into chunks with {str(contigs_per_file)} contigs per file. Writing files to {output_dirname} \n")
    records_chunk=[]
    for n, record in enumerate(screed.open(fasta)):
    #for record in screed.open(fasta):
        record.name = record.name.rsplit("\t", 1)[0]
        records_chunk.append(record)
        if n > 0 and n % contigs_per_file == 0:
            out_log.write(f"writing {str(n)}th contig to {str(filenum)} fasta file\n")
            outfile = os.path.join(output_dirname, prefix + "_" + str(filenum) + ".fa")
            with open(outfile, "w") as out:
                for record in records_chunk:
                    out.write(f">{record.name}\n{record.sequence}\n")
            records_chunk=[]
            filenum+=1
    #out_log.write(f"{str(i)} contigs written as individual fasta files\n")
    out_log.write(f"{str(n)} contigs written as {str(filenum)} individual fasta files\n")
