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
fasta_dir = config['fasta_dir']
num_chunks = 20000

# rather than using a checkpoint, let's:
# 1. Given orthogroup info and a folder with the split fastas: 
# 2. Split into num_chunks, which is the number of sketch jobs to run
# 3. write filelist for each chunk
# 4. sketch each chunk with --from-file, scaled=1
# 5. zip sketches from each chunk (by alpha ksize)
# 6. zip sketches from all chunks --> single zipfile (by alpha ksize)
# 7. sbts from the zip

## scaled = 1 for initial sketches and zipfiles
# protein k=7. k=10, dayhoff k=16
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
og_chunks={}

for i in range(0, num_chunks):
    og_chunks[i] = ORTHOGROUPS[i::num_chunks]

rule all:
    input:
        expand(f"{out_dir}/index/{basename}.{{ak}}-scaled{{scaled}}.sbt.zip", ak=alpha_ksizes, scaled=scaled)


localrules: write_chunked_fasta_filelists
rule write_chunked_fasta_filelists:
    output: expand(os.path.join(out_dir, 'og_chunks', f"{basename}.{{chunk_idx}}.txt"), chunk_idx = range(0, num_chunks))
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        runtime=30
    threads: 1
    run:
        for i in range(0, num_chunks):
            this_outfile = os.path.join(out_dir, "og_chunks", f"{basename}.{i}.txt")
            this_chunk = og_chunks[i] #ORTHOGROUPS[i::num_chunks]
            print(f"{i} chunk length:",len(this_chunk))
            with open(this_outfile, 'w') as out:
                for og in this_chunk:
                    fasta_path = os.path.join(fasta_dir, f"{basename}.{og}.fa")
                    if os.path.exists(fasta_path):
                        out.write(fasta_path + "\n")
                    else:
                        print(f"Missing fastafile! This should not happen {fasta_path}\n")



rule sourmash_sketch_chunks:
    input: 
        os.path.join(out_dir, 'og_chunks', f"{basename}.{{chunk_idx}}.txt"),
    output: 
        #expand(os.path.join(out_dir, "sigs", f"{basename}.{orthogroup}.fa.sig"), orthogroup=og[w.chunk_idx])
        os.path.join(out_dir, "og_chunk_sigs", f"{basename}.{{chunk_idx}}.sketched.txt"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000,#3000,
        runtime=6000,
    log: os.path.join(logs_dir, 'sketch_chunks', f"{basename}.{{chunk_idx}}.log")
    benchmark: os.path.join(logs_dir, 'sketch_chunks', f"{basename}.{{chunk_idx}}.benchmark")
    params:
        param_string=all_p_str,
        outdir=os.path.join(out_dir, "sigs"),
    conda:
        "conf/env/sourmash.yml" 
    #threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        sourmash sketch protein {params.param_string} --from-file {input} --outdir {params.outdir} 2> {log}
        touch {output}
        """

localrules: write_sketchlist
rule write_sketchlist:
    input: 
        #expand(os.path.join(out_dir, "og_chunk_sigs", f"{basename}.{{chunk_idx}}.sig"), chunk_idx= range(0,num_chunks))
        expand(os.path.join(out_dir, "og_chunk_sigs", f"{basename}.{{chunk_idx}}.sketched.txt"), chunk_idx = range(0, num_chunks)),
        #expand(os.path.join(out_dir, "sigs", "{orthogroup}.sig"), orthogroup=ORTHOGROUPS)
    output: 
        os.path.join(out_dir, f"{basename}.siglist.txt")
    resources: 
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        runtime=30
    threads: 1
    run:
        with open(str(output), 'w') as out:
            for og in ORTHOGROUPS:
                sketch_path = os.path.join(out_dir, "sigs", "{og}.sig")
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
