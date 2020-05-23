 

## eggnog

 eggnog 5.0: http://eggnog5.embl.de/download/eggnog_5.0/
 eggnog proteomes: http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa

 annotations tsv: http://eggnog5.embl.de/download/eggnog_5.0/e5.og_annotations.tsv
 tax id: http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv


## orthodb: https://www.orthodb.org/?page=filelist
    
`zgrep 1000373_1 orthodb/*` yields:
    ```
    odb10v1_all_fasta.tab.gz:>1000373_1:000000	1000373_1
    odb10v1_all_fasta.tab.gz:>1000373_1:000001	1000373_1
    odb10v1_all_og_fasta.tab.gz:>1000373_1:000001	1000373_1
    odb10v1_genes.tab.gz:1000373_1:000000	1000373_1	YP_005097973.1	RnQV1s4_gp1			11604942	structural protein
    odb10v1_genes.tab.gz:1000373_1:000001	1000373_1	RdRP_4@614:1065@YP_005097975.1	RdRP_4@614:1065@YP_005097975.1				RNA-directed RNA polymerase, luteovirus
    odb10v1_gene_xrefs.tab.gz:1000373_1:000000	RnQV1s4_gp1	NCBIgenename
    odb10v1_gene_xrefs.tab.gz:1000373_1:000001	IPR001795	InterPro
    odb10v1_level2species.tab.gz:10239	1000373_1	2	{10239,35325,1000373}
    odb10v1_OG2genes.tab.gz:84at35325	1000373_1:000001
    odb10v1_OG2genes.tab.gz:9167at10239	1000373_1:000001
    odb10v1_species.tab.gz:1000373	1000373_1	Rosellinia necatrix quadrivirus 1	GCF_000895895.1	1	2	C
    ```

So, I think we want to keep full `1000373_1:000000` name as the signature name in sourmash sbt. 
To do this, approximately:

    ```
    for record in screed.open("odb10v1_all_fasta.tab.gz"):
        signame = (record.name).rsplit("\t", 1)[0]
        mh = determine_appropriate_fresh_minhash(ksize, scaled, abundance, alphabet)
        mh.add_protein(record.sequence)
        sig = sourmash.Signature(mh, name=signame)
        leaf = SigLeaf(md5, sig)
        sbt.add_node(leaf)
    ```

Approach:
 
    1. grow sbts (but need to wait for resolution of https://github.com/dib-lab/sourmash/pull/994)
      - scaled=1, one signature per fasta record in odb10v1_all_fasta.tab.gz
    
    2. build test datasets:
        - 1. ~ positive control: exact fastas from orthodb odb10v1_all_fasta.tab.gz (should match 100%)
        - 2. less exact true positives: find same genes from other organisms. QfO peptides for now!
              - Find QfO peptides orthodb annotation, if it exists. if not: generate via orthodb tools.
        - 3. true negatives: use QfO ncRNA
              - try annotating these using orthodb methods, so we can compare where they get it wrong too (or, less likely, where there might be potential cds seqs within ncRNA)
    
    3. search/gather sbt with some genes we *should* be able to find (true positives)
        - how to present data? match to odb10v1_genes.tab.gz; return all info in csv (maybe gff3 later)
    
    3. Build variant of charcoal concept (~ gather using all query seqs --> build small LCA dbs --> per-contig identification)
    
    4. Compare to NCBI mapping approach used in nf-predictorthologs (sourmash search to, e.g. refseq sequences (any db allowed):
      `https://github.com/czbiohub/test-datasets/raw/predictorthologs/reference/ncbi_refseq_vertebrate_mammalian_ptprc_plus__np_only.fasta`
   
        ```
        sourmash index --ksize ${sourmash_ksize} --${sourmash_molecule} ${reference_proteome_sig.simpleName} ${reference_proteome_sig}
        ```
        Then search
        ```   
        sourmash search \\
            --containment \\
            --threshold 1e-100 \\
            --output ${csv_output} \\
            --ksize ${sourmash_ksize} \\
            --${sourmash_molecule} \\
            ${query_sig} \\
            ${reference_sbt_json}
        ```
        ** what is happening with the output? **

   5. Could my approach improve NCBI db use? Not sure we could directly build this sort of db with ncbi. Might end up with lots of unhelpful duplicates / it might just not be feasible. In order to work, this approach needs an sbt with gene-level (orthogroup-level) annotation mappings.
          - create NCBI databases in protein space, with lineage mappings?
          - think about whether it make sense to apply same concept
