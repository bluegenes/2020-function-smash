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

downloadable_accessions = [g for g in accessions if g.startswith("GC")] # 15124 with genome accessions
# rows with no genbank accession. Most have no entry in this column; Triticum aestivum has broken link
# ncbi has info for some of these. for each tax id, can run:
# scripts/ncbi-datasets assembly-descriptors tax-id 4565 --limit ALL to see info, including accession for representative genome
# Can parse the json output with jq to get just the accession(s):
# scripts/ncbi-datasets assembly-descriptors tax-id 5656 --limit ALL | jq '.datasets[].assembly_accession' -r
no_gc_accession = odb_species[~odb_species[ "genome-accession"].astype(str).str.startswith('GC', na=False)] #123 weird ones
#taxid2accession = no_gc_accession.set_index('ncbi-tax-id').to_dict()['genome-accession']
taxid2sciname = no_gc_accession.set_index('ncbi-tax-id').to_dict()['scientific-name']
really_no_genome_list = []

#def try_to_lookup_accession(taxid, sciname, refseq=False):
## this lookup works, but only prints accession to command line, does not store to a variable we can use.
## how to get around? Could print to a file. Otherwise...?
#    try:
#        shell("scripts/ncbi-datasets assembly-descriptors tax-id {taxid} --refseq | jq '.datasets[].assembly_accession' -r")  > accession
#        return accession
#    except:
#        pass
#    try:
#        shell("scripts/ncbi-datasets assembly-descriptors tax-id {taxid} --limit ALL | jq '.datasets[].assembly_accession' -r") >  accession
#        return accession
#    except:
#        pass
#    try: 
#        return shell("scripts/ncbi-datasets assembly-descriptors tax-name \"{sciname}\" --limit ALL | jq '.datasets[].assembly_accession' -r")
#    except:
#        return None

#for taxid, acc in taxid2accession.items():
#    accession= try_to_lookup_accession(taxid, taxid2sciname[taxid])
#    import pdb;pdb.set_trace()
#    if accession:
#        taxid2accession[taxid] = accession
#    else:
#        really_no_genome_list.append(taxid)





# write out tax info with headers :) 
odb_species.to_csv("odb10v1_species.all.csv", index=False)
no_gc_accession.to_csv("odb10v1_species.no_ncbi_download.csv", index=False)

rule all:
    input: 
        expand("orthodb/species_genomes/{genome}.fna.gz", genome = downloadable_accessions),
        expand("orthodb/taxid_genomes/{taxid}.{ext}", taxid = taxid2sciname.keys(), ext=["zip", "info.json", "accessions.txt"])


#alternate tool (cheers tereiter):
rule download_ncbi_datasets_tool:
    output: "scripts/ncbi-datasets"
    shell:
        """
        # linux version 
        wget -O {output} wget -O scripts/ncbi-datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets 
        chmod +x {output}
        """

# this rule will fail for species without any genomes available at ncbi
rule download_all_taxid:
    input: "scripts/ncbi-datasets"
    output: 
        info="orthodb/taxid_genomes/{taxid}.info.json",
        accession_info="orthodb/taxid_genomes/{taxid}.accessions.txt",
        zipfile="orthodb/taxid_genomes/{taxid}.zip"
    threads: 1
    resources:
        mem_mb= 2000,
        runtime= 60
    shell:
        """
        scripts/ncbi-datasets assembly-descriptors tax-id {wildcards.taxid} --limit ALL > {output.info}
        scripts/ncbi-datasets assembly-descriptors tax-id {wildcards.taxid} --limit ALL | jq '.datasets[].assembly_accession' -r > {output.accession_info}
        scripts/ncbi-datasets download assembly tax-id {wildcards.taxid} -f {output.zipfile}
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



