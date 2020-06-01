"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s grow-singleton-sbt.snakefile --configfile config/orthodb.yml --use-conda
"""

import os
import glob
import pandas as pd

out_dir = config.get("out_dir", "singleton_sbts")
logs_dir = os.path.join(out_dir, "logs")

dbInfo=config["databases"]

sbt_targets=[]

for database, info in dbInfo.items():    
    # build sbt targets
    for alpha, alphainfo in info["alphabet"].items():
        sbt_targets+=expand(os.path.join(out_dir,"{db}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"), db=database, alphabet=alpha, scaled=alphainfo["scaled"], k=alphainfo["ksizes"])

rule all:
    input: sbt_targets

rule grow_singleton_sbt:
    input:
        lambda w: dbInfo[w.database]["input_files"]
    output: 
        sbt=os.path.join(out_dir,"{database}.{alphabet}_scaled{scaled}_k{k}.sbt.zip"),
    threads: 1
    params:
        alpha= lambda w: w.alphabet.rsplit("translate_")[1] if w.alphabet.startswith("translate") else w.alphabet, # remove translate
        translate = lambda w: " --translate " if w.alphabet.startswith("translate") else "",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=60000,
    log: os.path.join(logs_dir, "grow-sbt", "{database}.{alphabet}_scaled{scaled}_k{k}.grow-singleton.log")
    benchmark: os.path.join(logs_dir, "grow-sbt", "{database}.{alphabet}_scaled{scaled}_k{k}.grow-singleton.benchmark")
    conda: "envs/sbt-env.yml"
    shell:
        # escaped quotes allows for "*" in input dir
        #python scripts/grow-sbtmh.py \"{params.input_dir}\" --input-is-directory --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {params.alpha} --subset-csv {params.subset_csv} --subset-info-colname {params.subset_info_colname} {params.translate} 2> {log}
        """
        python scripts/grow-sbtmh.py {input} --sbt {output.sbt} --ksize {wildcards.k} --scaled {wildcards.scaled} --alphabet {params.alpha} --singleton {params.translate} 2> {log}
        """
