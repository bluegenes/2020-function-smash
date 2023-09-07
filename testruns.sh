


#python scripts/grow-sbtmh.py test_grow_sbt --sbt test_grow_sbt_subset/sbt-7paths_dayhoff_k13.sbt.zip --ksize 13 --input-is-directory --scaled 1 --alphabet dayhoff --info-colname accession --subset-csv gtdb_evol_groups_first_seven.tsv

python grow-sbtmh.py orthodb/odb10v1_all_fasta.tab.gz --sbt function_sbts/odb10v1_all_fasta.dayhoff_k13_scaled1.sbt.zip --ksize 13 --scaled 1 --alphabet dayhoff
