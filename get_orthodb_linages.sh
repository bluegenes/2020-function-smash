
# test small
#python scripts/get_lineage_from_taxid.py test_data/odb10v1_species.head200.csv --output-csv test_lin.csv --charcoal-csv test_lin.charcoal.csv


# full
python scripts/get_lineage_from_taxid.py odb10v1_species.all.csv --output-csv odb10v1_species.all.ncbi-lineage.tsv --charcoal-csv odb10v1_species.all.charcoal.csv


