 #!/bin/bash

#snakemake -s download-ncbi.snakefile  --profile farm --cluster-config cluster_config.yml --jobs 30

#snakemake -s download-ncbi.snakefile  --profile default  --jobs 7 -n

#snakemake -s grow-singleton-sbt.snakefile --profile default --configfile config/orthodb.yml

#snakemake -s grow-singleton-sbt.snakefile --profile farm  --cluster-config cluster_config.yml --configfile config/orthodb.yml --jobs 1

#snakemake -s large-fasta-sigs-to-singleton-sbt.snakefile --profile farm --cluster-config cluster_config.yml --configfile config/orthodb.yml --jobs 1

#snakemake -s large-fasta-sigs-to-singleton-sbt.snakefile --profile default  --configfile config/orthodb.yml --jobs 5


snakemake -s index-by-orthogroup.snakefile --profile farm --cluster-config cluster_config.yml --configfile config/orthodb.yml --jobs 30

#snakemake -s index-by-orthogroup.snakefile --configfile config/orthodb.yml --profile default --cluster-config cluster_config.yml --jobs 1 --restart-times 0



