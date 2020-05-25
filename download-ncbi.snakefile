"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  download-ncbi.snakefile -k 
"""

import os
import re
import pandas as pd

#####
#odb10v1_species.tab
#1.  NCBI tax id
#2.  Ortho DB individual organism id, based on NCBI tax id
#3.  scientific name inherited from the most relevant NCBI tax id
#4.  genome asssembly id, when available
#5.  total count of clustered genes in this species
#6.  total count of the OGs it participates
#7.  mapping type, clustered(C) or mapped(M)

header=["ncbi-tax-id", "orthodb-organism-id", "scientific-name", "genome-accession", "num-clustered-genes", "num-orthogroups", "clustered-or-mapped"]
odb_species = pd.read_csv("orthodb/odb10v1_species.tab.gz", sep="\t", names=header)
accessions=odb_species["genome-accession"].dropna().astype(str).tolist()

downloadable_accessions = [g for g in accessions if g.startswith("GC")][:1] # 15124 with genome accessions
# rows with no genbank accession. Most have no entry in this column; Triticum aestivum has broken link
no_gc_accession = odb_species[~odb_species[ "genome-accession"].astype(str).str.startswith('GC', na=False)] #123 weird ones
odb_species.to_csv("odb10v1_species.all.csv", index=False)
no_gc_accession.to_csv("odb10v1_species.no_ncbi_download.csv", index=False)


rule all:
    input: expand("orthodb/species_genomes/{genome}.fna.gz", genome = downloadable_accessions)

#alternate tool (cheers tereiter):
rule download_ncbi_datasets_tool:
    output: "scripts/ncbi-datasets"
    shell:
        """
        # linux version 
        wget -O {output} wget -O scripts/ncbi-datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets 
        chmod +x {output}
        """


rule ncbi_datasets_download:
    input: "scripts/ncbi-datasets"
    output: "orthodb/species_genomes/{genome}.zip"
    threads: 1
    resources:
        mem_mb= 2000,
        runtime= 60
    shell:
        """
        scripts/ncbi-datasets download assembly {wildcards.genome} -f {output}
        """


rule unzip_genomes:
    input: "orthodb/species_genomes/{genome}.zip"
    output: "orthodb/species_genomes/{genome}.fna.gz"
    resources:
        mem_mb= 3000,
        runtime= 60
    shell:
        """
        unzip -p {input} ncbi_dataset/data/{wildcards.genome}/*.fna | gzip -9 > {output}
        """

# sigh, why you no work, ncbi-genome-downloader?

# test genome GCF_001500715.1
#rule download_by_accession:
#    input: "scripts/datasets"
#    output: "orthodb/species_genomes/{genome}.fna.gz"
#    params:
#        out_folder="orthodb/species_genomes"
#    conda: ncbi-genome-downloader.yml
#    shell:
#        """
#        ncbi-genome-download all --assembly-accessions {genome} --output-folder {params.out_folder} --flat-output --format 'fasta'
#        """



