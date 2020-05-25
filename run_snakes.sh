 #!/bin/bash

#snakemake -s download-ncbi.snakefile  --profile farm --cluster-config cluster_config.yml --jobs 5

snakemake -s download-ncbi.snakefile  --profile default  --jobs 7


