from pathlib import Path
import subprocess
from node import Node

# Inputs (assumed, will change if needed):
#   - list[str] containing only the sequences, no labels
#   - str with model type (eg. "LG")

fasta_path = Path("homologous_bglb_sequences.fasta")
output_dir = Path("iqtree_output")

def create_fasta(sequences: list[str]) -> None:
    with fasta_path.open("w") as f:
        for i in range(len(sequences)):
            f.write(f">protein_{i}\n")
            f.write(f"{sequences[i]}\n")

def run_iqtree(model: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / "bglb_asr"

    command = ["iqtree2", "-s", fasta_path, "-st", "AA", "-m", model, "-asr",
        "-pre", str(prefix)
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)

def extract_sequences() -> dict[str, str]:
    state_file = output_dir / "bglb_asr.state"
    sequences = {}
    with state_file.open("r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            # Fields[2] should have the state
            if fields[2] == "State": # This is the header line
                continue
            if fields[0] not in sequences.keys():
                sequences[fields[0]] = ""
            sequences[fields[0]] += fields[2]
    return sequences

def build_tree(node_dict: dict[str, str], protein_seqs: list[str]) -> list[Node]
    tree_file = output_dir / "bglb_asr.treefile"
    with tree_file.open("r") as f:
        newick = next(f)

    i = 0
    stack = []
    tree = {}
    while i < len(newick):
        if newick[i] == "(":
            name = newick[i+1, newick.find(":", i)]
            stack.append(name)
            if "protein" in name:
                num = int(name.split("_")[1])
                seq = protein_seqs[num]
            else:
                seq = node_dict[name]
            tree[name] = Node(seq, None)
            i = newick.find(":", i) + 1

        elif newick[i] == ")":
            # the last two nodes appended will be the children of the root
            root_name = newick[i+1, newick.find(":", i)]
            child_1 = stack.pop()
            child_2 = stack.pop()

if __name__ == "__main__":
    test_hom = ["MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNTD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNTE",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVSTD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVSTE",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNAD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNAE"
    ]
    #create_fasta(test_hom)
    #run_iqtree("LG")
    seqs = extract_sequences()
    print(seqs)
