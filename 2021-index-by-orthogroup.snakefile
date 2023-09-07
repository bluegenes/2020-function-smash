"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s 2021-index-by-orthogroup.snakefile
"""

import os
import glob
import pandas as pd


configfile: "conf/orthodb.yml"
out_dir = config.get("output_dir", "outputs")
logs_dir = os.path.join(out_dir, "logs")
basename = config['basename']

# rather than using a checkpoint, let's:
# 1. read through orthogroup tab file, gget all orthogroups
# 2. split the fasta, expecting one file per orthogroup
# 3. sketch each orthogroup, scaled=1
# 4. write_sketchlist
# 5. zip_sketches
# 6. sbts from the zips

## scaled = 1 for initial zketches and zipfiles
# protein k=7. k=10, dayhoff k=16
# temporary split fastas? 
alphabet_info = config.get("alphabet_info", {'protein': [7], 'dayhoff': [16]})
scaled=config.get("scaled", [1])
minscaled = min(scaled)

p_str = []
alpha_ksizes=[]

for alpha, ksizes in alphabet_info.items():
    ak=expand(f"{alpha}-k{{k}}", k=ksizes)
    alpha_ksizes += ak
    # sourmash param strings
    kinf = expand("k={k}", k=ksizes)
    kstr = ",".join(kinf)
    this_pstr = f"-p {alpha},{kstr},scaled={minscaled},abund"
    p_str.append(this_pstr)

all_p_str = " ".join(p_str)
print(all_p_str)
print(alpha_ksizes)


# get all orthogroups
ORTHOGROUPS=[]
og_info = config["og_info"]
ORTHOGROUPS = pd.read_csv(og_info, usecols=[0], sep='\t', header=None, squeeze=True).tolist() # squeeze forces --> series instead of DF, so can use tolist()

rule all:
    input:
        expand(f"{out_dir}/index/{basename}.{{ak}}-scaled{{scaled}}.sbt.zip", ak=alpha_ksizes, scaled=scaled)

checkpoint split_by_orthogroup:
    input: 
        fasta= config["input_fasta"],
        og2genes= config["gene_mapping"]
    output: 
        #temp(expand(os.path.join(out_dir, 'split_fasta', f"{basename}.{{orthogroup}}.fa"), orthogroup = ORTHOGROUPS))
        directory(os.path.join(out_dir, 'split_fasta')),
    params:
        prefix= basename,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    log: os.path.join(logs_dir, "split_by_orthogroup", f"{basename}_split.log")
    benchmark: os.path.join(logs_dir, "split_by_orthogroup", f"{basename}_split.benchmark")
    conda: "conf/env/sourmash.yml" 
    script: os.path.join("scripts", "split-by-orthogroup-dict.py")


rule sourmash_sketch:
    input: os.path.join(out_dir, 'split_fasta', f"{basename}.{{orthogroup}}.fa"),
    output: 
        #temp(os.path.join(out_dir, "sigs", "{orthogroup}.sig"))
        os.path.join(out_dir, "sigs", "{orthogroup}.sig")
    params:
        param_string =all_p_str
    conda: "conf/env/sourmash.yml" 
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        runtime=60,
    threads: 1
    shell:
        """
        sourmash sketch protein {params.param_string} -o {output} --name {wildcards.orthogroup} {input}
        """

localrules: write_sketchlist

# make the checkpoint work
def checkpoint_output_split_by_orthogroup(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.split_by_orthogroup.get(**wildcards).output[0]    
    file_names = expand(f"{out_dir}/sigs/{{orthogroup}}.sig", orthogroup = ORTHOGROUPS)  # use explicit orthogroups rather than glob_wildcards 
                        #orthogroup = glob_wildcards(os.path.join(checkpoint_output, "{base}.fa")).pfam)
    return file_names


localrules: write_sketchlist
rule write_sketchlist:
    input: 
        #expand(os.path.join(out_dir, "sigs", "{orthogroup}.sig"), orthogroup=ORTHOGROUPS)
        checkpoint_output_split_by_orthogroup
    output: 
        os.path.join(out_dir, f"{basename}.siglist.txt")
    resources: 
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        runtime=30
    threads: 1
    run:
        with open(str(output), 'w') as out:
            for inF in input:
                sketch_path = str(inF)
                if os.path.exists(sketch_path):
                    out.write(sketch_path + "\n")
                else:
                    print(f"Missing sketchfile! This should not happen {sketch_path}\n")


rule zip_sketches:
    input:
        siglist=os.path.join(out_dir, f"{basename}.siglist.txt"),
        #sigs=expand(os.path.join(out_dir, "sigs", "{orthogroup}.sig"), orthogroup=ORTHOGROUPS)
        #sigs=checkpoint_output_split_by_orthogroup # don't need this here if sigs not temporary
    output: 
        os.path.join(out_dir, "zip", f"{basename}.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip")
    log: 
        os.path.join(logs_dir, "zip", f"{basename}.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip.log")
    params:
        alpha_cmd= lambda w: f"--{w.alphabet}",
    conda: "conf/env/sourmash.yml" 
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5000,
        runtime=600,
    threads: 1
    shell:
        '''
        sourmash sig cat --from-file {input} {params.alpha_cmd} --ksize {wildcards.ksize} -o {output} 2> {log}
        '''


rule sbt_index:
    input:
         os.path.join(out_dir, "zip", f"{basename}.{{alphabet}}-k{{ksize}}-scaled{minscaled}.zip")
    output: 
         os.path.join(out_dir, "index", f"{basename}.{{alphabet}}-k{{ksize}}-scaled{{scaled}}.sbt.zip")
    params:
        alpha_cmd= lambda w: f"--{w.alphabet}",
    log: 
        os.path.join(logs_dir, "index", f"{basename}.{{alphabet}}-k{{ksize}}-scaled{{scaled}}.sbt.log")
    conda: "conf/env/sourmash.yml" 
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 50000,
        runtime=600,
    threads: 1
    shell:
        '''
        sourmash index {output} {input} {params.alpha_cmd} --ksize {wildcards.ksize} --scaled {wildcards.scaled} 2> {log}
        '''
