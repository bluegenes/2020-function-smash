"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s singleton-sig-sbt.snakefile --use-conda
"""

import os
import glob
import pandas as pd

out_dir = config.get("out_dir", "output-singleton-smash")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"
compute_dir = os.path.join(out_dir, "compute")
index_dir = os.path.join(out_dir, "index")

index_targets = []

sampleInfo=config["databases"]

for sample, info in sampleInfo.items():    
    for alpha, alphainfo in info["alphabet"].items():
        index_targets+=expand(os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.singleton.sbt.zip"), sample=sample, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: index_targets

# chunk large fasta to parallelize signature generation!
localrules: split_fasta

checkpoint split_fasta:
    input: lambda w: sampleInfo[w.sample]["input_fasta"]
    output:
        directory(os.path.join(out_dir, "{sample}_fastasplit")),
    params:
        prefix=lambda w: w.sample,
        contigs_per_file=1000, # ~40 million reads in orthodb --> will create 40k files
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "split_fasta", "{sample}_split.log")
    benchmark: os.path.join(logs_dir, "split_fasta", "{sample}_split.benchmark")
    conda: os.path.join(envs_dir, "khmer-env.yml") # requires screed
    script: os.path.join("scripts", "split-fasta.py")

# for protein signatures, multipy by 3 if necessary before calculating signature (sourmash v3.x)
ksize_multiplier = {"nucleotide": 1, "protein": 3, "dayhoff": 3, "hp":3, "translate_protein": 3, "translate_dayhoff": 3, "translate_hp": 3}
# "rna" is not sourmash cli friendly
moltype_map = {"nucleotide": "dna", "protein":"protein", "dayhoff": "dayhoff", "hp":"hp", "translate_protein": "protein", "translate_dayhoff": "dayhoff", "translate_hp": "hp"}

# compute sigs to general compute directory, regardles of the sbt sample name (so can reuse sigs for other sbts, if desired)
rule sourmash_compute_dna:
    input: os.path.join(out_dir, "{sample}_fastasplit", "{sample}_{tnum}.fa") 
    output: os.path.join(compute_dir, "{sample}.dna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{sample}.dna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}.dna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

rule sourmash_compute_protein:
    input: os.path.join(out_dir, "{sample}_fastasplit", "{sample}_{tnum}.fa") 
    output: os.path.join(compute_dir, "{sample}.protein", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=True,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{sample}.protein", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.protein.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}.protein", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.protein.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"

rule sourmash_compute_rna:
    input: os.path.join(out_dir, "{sample}_fastasplit", "{sample}_{tnum}.fa") 
    output: os.path.join(compute_dir, "{sample}.rna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.sig")
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: moltype_map[w.alphabet],
        input_is_protein=False,
        track_abundance=True,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{sample}.rna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.rna.compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}.rna", "{sample}_{tnum}_{alphabet}_scaled{scaled}_k{k}.rna.compute.benchmark")
    conda: "envs/sourmash3.3.yml"
    script: "scripts/sourmash-compute.wrapper.py"


def aggregate_sigs(w):
    """
    aggregate the file names of the unknown number of files generated at the split_fasta step
    """
    input_type = sampleInfo[w.sample]["input_type"]
    assert input_type in ["protein", "dna", "rna"]
    sig_dir = os.path.join(compute_dir, f"{w.sample}.{input_type}")
    # get checkpoint output from split_fasta
    checkpoint_output=checkpoints.split_fasta.get(sample=w.sample).output[0]
    # build pattern we'll be using below. globbing these files will populate tnum wildcard
    fasta_pattern = os.path.join(checkpoint_output, f"{w.sample}_{{tnum}}.fa")
    # build signature path, including all fixed variables (sample, alphabet, scaled, ksize)
    sigpath= os.path.join(sig_dir, f"{w.sample}_{{tnum}}_{w.alphabet}_scaled{w.scaled}_k{w.k}.sig")
    # expand sigpath to sigfiles by globbing for tnum wildcard in the checkpoint output (fasta pattern)
    sigfiles=expand(sigpath, tnum=glob_wildcards(fasta_pattern).tnum)
    #print(sigfiles)
    return sigfiles 


rule grow_singleton_sbt:
    input: aggregate_sigs,
    output:
        sbt=os.path.join(index_dir,"{sample}.{alphabet}_scaled{scaled}_k{k}.singleton.sbt.zip"),
    threads: 1
    params:
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
        input_type = lambda w: sampleInfo[w.sample]["input_type"],
        ksize = lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "grow-singleton-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.log")
    benchmark: os.path.join(logs_dir, "grow-singleton-sbt", "{sample}.{alphabet}_scaled{scaled}_k{k}.grow.benchmark")
    conda: "envs/sbt-env.yml"
    shell:
        """
        python scripts/grow-sbtmh.py {compute_dir}/{params.input_type}/*_{params.alpha}_scaled{wildcards.scaled}_k{wildcards.k}.sig \
        --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {params.alpha} --singleton {params.translate} 2> {log}
        """
