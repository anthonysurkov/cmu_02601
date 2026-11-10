
# (3a) IQ-TREE ModelFinder wrapper

# Input:  list of homologous BglB sequences (strings)
# Output: the name of the best amino-acid substitution model, e.g. "LG+F+I+G4"

# had to install iqtree2 and mafft  

import subprocess
import os


def find_model(sequences, workdir="modelfinder_out"):
    """Runs ModelFinder on the given sequences and returns the best model name."""

    # make a folder to keep all the IQ-TREE files in
    os.makedirs(workdir, exist_ok=True)
    raw_file = os.path.join(workdir, "homologs.fasta")
    aln_file = os.path.join(workdir, "aligned.fasta")
    prefix = os.path.join(workdir, "run")

    #  writing the sequences out to a FASTA file
    with open(raw_file, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(">seq" + str(i) + "\n")
            f.write(seq.strip().upper() + "\n")


    #blast results are different lengths 
    #need to align, using mafft since modelfinder needs all same length 
    with open(aln_file, "w") as out:
        subprocess.run(["mafft", "--auto", raw_file], stdout=out, check=True)

    #  run IQ-TREE in ModelFinder-only mode -m MF means just pick a model, no tree 
    subprocess.run([
        "iqtree",
        "-s", aln_file,         # the alignment
        "-m", "MF",             # ModelFinder only
        "--seqtype", "AA",       # amino acids not DNA
        "--prefix", prefix,     # name all the output files "run.*"
        "--redo",                # overwrite old results if we run it again
    ], check=True)

    # 4. read the report file and find the line with the winning model on it.
    #    The line looks like: "Best-fit model according to BIC: LG+F+I+G4"
    with open(prefix + ".iqtree") as f:
        for line in f:
            if "Best-fit model according to" in line:
                model = line.split(":")[1].strip()
                return model

    # if we get here something went wrong
    raise RuntimeError("Couldn't find the best-fit model. Check " + prefix + ".log")


# testing 
if __name__ == "__main__":
    test_seqs = [
        "MKFPKDFMWGTATAAYQIEGA",
        "MKFPEDFLWGTATAAYQIEGA",
        "MRFPKDFMWGAATAAYQVEGA",
    ]
    print(find_model(test_seqs))

    #type hint 
    #two lines between classes and imports 
    #self.directory and paths.py data dir 
