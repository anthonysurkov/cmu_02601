#Module 2: BLAST API wrapper.
#Takes the canonical BglB sequence and returns homologous sequences

from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from .paths import DATA_DIR


BLAST_PROGRAM = "blastp" #protein against protein 
BLAST_DATABASE = "swissprot"  #chosen dataset which is one half of uniprot 
MAX_HITS = 250  # blast query over database for BgIB homologs retreives 250 hits 
E_VALUE_CUTOFF = 1e-10  # strict setting to only get the good matches 
MIN_IDENTITY = 0.35 # a hit must match your BglB at 35% of aligned positions 
CACHE_DIR = DATA_DIR / "blast"


#creating filename with residues and seq length 
def _get_cache_path(sequence: str) -> Path: 
    tag = f"{len(sequence)}_{sequence[:10]}"
    return CACHE_DIR / f"blast_{tag}.xml" 


#getting blast results from ncbi and saving them to disk 
def _run_blast(sequence: str, cache_path: Path) -> None:
    result_handle = NCBIWWW.qblast(
        program=BLAST_PROGRAM,
        database=BLAST_DATABASE,
        sequence=sequence,
        hitlist_size=MAX_HITS,
        expect=E_VALUE_CUTOFF,
    )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as out_file:
        out_file.write(result_handle.read())
    result_handle.close()


#getting the percent identity for one alignment 
# hsp is blast object that shows one aligned region between query and database (High scoring Segment Pair)
#hsp.identities is num positions in alignment holding same amino acid in both sequences
def _get_identity(hsp, alignment_length: int) -> float:
    if alignment_length == 0:
        return 0.0
    return hsp.identities / alignment_length


#turning the saved xml file into list of sequences 
def _parse_blast(cache_path: Path) -> list[str]:
    sequences = []  # grows to at most MAX_HITS

    with open(cache_path) as xml_file: #loading the file 
        record = next(NCBIXML.parse(xml_file))

    for alignment in record.alignments:  # up to 250 alignments
        hsp = alignment.hsps[0]  # best scoring local alignment only
        identity = _get_identity(hsp, hsp.align_length)

        if hsp.expect > E_VALUE_CUTOFF or identity < MIN_IDENTITY: #filter 
            continue

        # gaps are an artifact of the pairwise alignment in module 4 
        # realigns everything, so store the plain hit sequence.
        sequences.append(hsp.sbjct.replace("-", ""))

    return sequences


#given the reference BglB sequence, returns that sequence plus its BLAST homologs
# querying NCBI only if not alr have the results cached
def get_homologs(sequence: str, use_cache: bool = True) -> list[str]: 
    cache_path = _get_cache_path(sequence)

    if not (use_cache and cache_path.exists()):
        _run_blast(sequence, cache_path)

    homologs = _parse_blast(cache_path) 
    print(f"kept {len(homologs)} hits out of {MAX_HITS} requested.")

    return [sequence] + homologs


#reading in a fasta file and getting the seq as a string 
def read_fasta(fasta_path: Path) -> str:
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return str(record.seq)


if __name__ == "__main__":
    canonical = read_fasta(DATA_DIR / "bglb_canonical.fasta") #reading reference file 
    print(f"Query length: {len(canonical)} residues.") #checking to make sure residues are around 400-500

    results = get_homologs(canonical)
    print(f"Total sequences for ASR: {len(results)}") #running full pipeline and getting final count 