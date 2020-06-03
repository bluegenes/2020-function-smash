import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
from sourmash.sbtmh import SigLeaf
import pandas as pd

def try_reading_csv(groups_file):
    # autodetect format
    samples=""
    if '.tsv' in groups_file or '.csv' in groups_file:
        separator = '\t'
        if '.csv' in groups_file:
            separator = ','
        try:
            samples = pd.read_csv(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in groups_file:
        try:
            samples = pd.read_excel(groups_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {groups_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples


def try_reading_fasta_file(fasta_file):
    records=[]
    try:
        records=screed.open(fasta_file)
    except:
        # log the error somehow
        sys.stderr.write(f"\ncannot read fasta file {fasta_file}\n")
    return records

def create_sbt_or_load_existing(tree_file, load_existing=False):
    # hmm.. adding and overwriting seems complicated but managed? see sourmash/sbt_storage.py
    if load_existing:
        try:
            sbt = sourmash.load_sbt_index(tree_file)
        except:
            sys.stderr.write(f"\ncannot load sbt file at {tree_file}\n")
    else:
        sbt = sourmash.create_sbt_index()
    return sbt

def determine_appropriate_fresh_minhash(alphabet, ksize, scaled_val, ignore_abundance=False):
    # default behavior is to track abundance
    abund = not ignore_abundance
    if alphabet == "nucleotide":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=False)
    elif alphabet == "protein":
        k=ksize*3 ## need to multiply bt 3 to get same ksize, bc add_protein method does k/3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=False)
    elif alphabet == "dayhoff":
        k=ksize*3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=True, hp=False)
    elif alphabet == "hp":
        k=ksize*3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=True)
    return mh

def load_or_generate_sig_from_file(input_file, alphabet, ksize, scaled, ignore_abundance=False, translate=False):
    sig=""
    if input_file.endswith(".sig"):
        # do I want to enable multiple sigs per file here?
        sig = sourmash.load_one_signature(input_file, ksize=ksize)
    else:
        # read file and add sigs
        records = try_reading_fasta_file(input_file)
        # build signature name from filename .. maybe just keep filename?
        #signame = os.path.basename(input_file.rsplit("_", 1)[0])
        # start with fresh minhash
        mh = determine_appropriate_fresh_minhash(alphabet, ksize, scaled, ignore_abundance)
        if records:
            for record in records:
                if alphabet == "nucleotide" or translate:
                    mh.add_sequence(record.sequence, force=True)
                else:
                    mh.add_protein(record.sequence)
            # minhash --> signature, using filename as signature name ..i think this happens automatically if don't provide name?
            sig = sourmash.SourmashSignature(mh, name=os.path.basename(input_file))
    return sig

def add_singleton_sigs(sbt, input_file, ksize, scaled, alphabet, ignore_abundance=False, translate=False):
    if input_file.endswith(".sig"):
        # maybe this? haven't tested!!
        sigs = sourmash.signature.load_signatures(input_file, ksize=ksize, select_moltype=alphabet)
        for sig in sigs:
            if sig.minhash:
                leaf = SigLeaf(sig.md5sum(), sig)
                sbt.add_node(leaf)
        # loop through and add each to sbt
    else:
        # read file and add sigs
        records = try_reading_fasta_file(input_file)
        # start with fresh minhash
        if records:
            for n, record in enumerate(records):
                signame = (record.name).rsplit("\t", 1)[0]
                if n % 10000 == 0:
                    sys.stderr.write(f"... building {n}th sig, {signame}\n")

                mh = determine_appropriate_fresh_minhash(alphabet, ksize, scaled, ignore_abundance)
                if alphabet == "nucleotide" or translate:
                    mh.add_sequence(record.sequence, force=True)
                else:
                    mh.add_protein(record.sequence)
            # minhash --> signature
                sig = sourmash.SourmashSignature(mh, name=signame)
                if sig.minhash:
                    leaf = SigLeaf(sig.md5sum(), sig)
                    sbt.add_node(leaf)
    return sbt


def collect_input_files_from_dir(input_dir, subset_csv=None, subset_info_colname="accession"):
    """
    Scan a directory for input fasta or signature files.
    If subset_csv is provided, glob the directory for the pattern provided in the
    subset_info_colname column.
    """
    input_files=[]
    if subset_csv:
        try:
            match_strings = set((try_reading_csv(subset_csv))[subset_info_colname].tolist())
            for m in match_strings:
                # to do: if have filename column, could directly use filenames rather than globbing
                input_files+=glob.glob(os.path.join(input_dir, f"*{m}*"))
        except:
            sys.stderr.write(f"can't collect input files from dir {input_dir} using subset_csv {subset_csv} column {subset_info_colname}")
            sys.exit()
    else:
        input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)]
    if not input_files:
        sys.stderr.write(f"can't collect input files from dir {input_dir}")
        sys.exit()
    return input_files


def grow_singleton_sbt(args):

    # find input files if necessary
    input_files = args.input_files
    if args.input_is_directory:
        input_files = collect_input_files_from_dir(input_files[0], args.subset_csv, args.subset_info_colname)

    # protein ksize acrobatics
    ksize = args.ksize
    alphabet = args.alphabet
    if alphabet == "protein" or alphabet == "dayhoff" or alphabet == "hp":
        ksize = ksize*3

    # load or create sbt
    sbt = create_sbt_or_load_existing(args.sbt, args.load_existing_sbt)

    # iterate through input files; add to sbt
    for n, filename in enumerate(input_files):
        # swipe some handy progress reporting code from titus:
        #if n % 100 == 0:
        sys.stderr.write(f"... loading {filename} file {n} of {len(input_files)}\n")
        sbt = add_singleton_sigs(sbt, filename, ksize, args.scaled, args.alphabet, args.ignore_abundance, args.translate)

    # save the tree
    sbt.save(args.sbt)

def grow_sbt(args):
    # find input files if necessary
    input_files = args.input_files
    if args.input_is_directory:
        input_files = collect_input_files_from_dir(input_files[0], args.subset_csv, args.subset_info_colname)
    ksize = args.ksize
    alphabet = args.alphabet
    if alphabet == "protein" or alphabet == "dayhoff" or alphabet == "hp":
        ksize = ksize*3
    # create or load sbt
    sbt = create_sbt_or_load_existing(args.sbt, args.load_existing_sbt)
    # iterate through input files; add to sbt
    for n, filename in enumerate(input_files):
        # swipe some handy progress reporting code from titus:
        if n % 100 == 0:
            sys.stderr.write(f"... loading {filename} file {n} of {len(input_files)}\n")

        # build or load signature from file

        sig = load_or_generate_sig_from_file(filename, alphabet, ksize, args.scaled, args.ignore_abundance, args.translate)
        # add to sbt
        if sig: # is this necessary?
            if sig.minhash:
                leaf = SigLeaf(sig.md5sum(), sig)
                sbt.add_node(leaf)

    # save the tree
    sbt.save(args.sbt)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("--sbt")
    p.add_argument("--ksize", type=int, default=31)
    p.add_argument("--scaled", type=int, default=1000)
    p.add_argument("--alphabet", default="nucleotide", help="desired alphabet for sourmash signatures. options: nucleotide, protein, dayhoff, or hp")
    p.add_argument("--translate", action="store_true", help="translate nucleotide input to the desired alphabet")
    p.add_argument("--input-is-directory", action="store_true", help="the input is a directory, rather than a file or series of files")
    p.add_argument("--subset-csv", default=None, help="provide a csv with info for choosing a subset of files from the input dir.")
    p.add_argument("--subset-info-colname", default="accession", help="specify the column name in the csv that will be used to glob files.")
    p.add_argument("--ignore-abundance", action="store_true", help="do not store abundance information for signatures")
    p.add_argument("--load-existing-sbt", action="store_true", help="load the sbt file if it exists. Turns off default overwrite of sbt file with new sbt.")
    p.add_argument("--singleton", action="store_true", help="with fasta inputs, add one signature per fasta entry, rather than one per file")
    args = p.parse_args()
    if args.singleton:
        sys.exit(grow_singleton_sbt(args))
    else:
        sys.exit(grow_sbt(args))
#to do: sourmash throws some warnings(?) when the sbt file exists already
