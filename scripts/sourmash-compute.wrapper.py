"""Snakemake wrapper for sourmash compute."""

__author__ = "Lisa K. Johnson and N. Tessa Pierce"
__copyright__ = "Copyright 2019, Lisa K. Johnson and N. Tessa Pierce"
__email__ = "ljcohen@ucdavis.edu and ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")
scaled = snakemake.params.get("scaled", "2000")
singleton = snakemake.params.get("singleton", False)
k = snakemake.params.get("k", "31")
compute_moltypes = snakemake.params.get("compute_moltypes", "dna")
input_is_protein = snakemake.params.get("input_is_protein", False)
track_abundance = snakemake.params.get("track_abundance", True)

k = [k] if (isinstance(k, str) or isinstance(k, int)) else k

compute_moltypes = [compute_moltypes] if isinstance(compute_moltypes, str) else compute_moltypes
moltype_cmd = ""
for moltype in compute_moltypes:
    if not moltype in ["dna", "protein", "dayhoff", "hp"]:
        raise ValueError(f"Can't compute signature for moltype {moltype}. Valid moltypes are 'dna', 'protein', 'dayhoff' or 'hp'")

#nucl_only_ksizes = []
prot_moltypes= ["protein", "dayhoff", "hp"]

if any((True for x in prot_moltypes if x in compute_moltypes)): #and not input_is_protein:
    for ksize in k:
        # if ksize not divisible by 3
        if (int(ksize) % 3 != 0):
            raise ValueError(f"some of your ksizes are not divisible by 3, can't translate to produce protein/dayhoff/hp sigs. If your input is protein, please set 'input_is_protein' to 'True'")

k = ",".join(map(str, k))

moltype_cmd = " --" + " --".join(compute_moltypes)
if input_is_protein:
    moltype_cmd += " --input-is-protein --no-dna "
abund_cmd = ""
if track_abundance:
    abund_cmd = " --track-abundance "
singleton_cmd=""
if singleton:
    singleton_cmd = " --singleton "


shell(
    "sourmash compute {moltype_cmd} {abund_cmd} {singleton_cmd} --scaled {scaled} -k {k} {snakemake.input} -o {snakemake.output} -p {snakemake.threads} {extra} {log}"
)
