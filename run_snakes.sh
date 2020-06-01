 #!/bin/bash

snakemake -s download-ncbi.snakefile  --profile farm --cluster-config cluster_config.yml --jobs 30

#snakemake -s download-ncbi.snakefile  --profile default  --jobs 7 -n

#snakemake -s grow-singleton-sbt.snakefile --profile default --configfile config/orthodb.yml

#snakemake -s grow-singleton-sbt.snakefile --profile farm  --cluster-config cluster_config.yml --configfile config/orthodb.yml --jobs 1

