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

def build_tree(node_dict: dict[str, str], protein_seqs: list[str]) -> list[Node]:
    tree_file = output_dir / "bglb_asr.treefile"
    with tree_file.open("r") as f:
        newick = next(f)

    i = 0
    stack = {}
    tree = {}
    stack_index = -1
    while i < len(newick):
        if newick[i] == "(":
            name = newick[i+1:newick.find(":", i)]
            stack_index += 1
            stack[stack_index] = [name]
            if "protein" in name:
                num = int(name.split("_")[1])
                seq = protein_seqs[num]
                tree[name] = seq
            else:
                seq = node_dict[name]
                tree[name] = Node(seq)
            i = newick.find(":", i) + 1

        elif newick[i] == ")":
            # the last two nodes appended will be the children of the root
            root_name = newick[i+1:newick.find(":", i)]
            root_name = root_name.strip(";")
            children = []
            for j in stack[stack_index]:
                children.append(tree[j])
            tree[root_name] = Node(node_dict[root_name], children)
            if ":" in newick[i+1:]:
                i = newick.find(":", i) + 1
                stack_index -= 1
                stack[stack_index].append(root_name)
            else: # Done
                i = len(newick)

        elif newick[i] == "," and newick[i+1].isalpha():
            name = newick[i+1:newick.find(":", i)]
            name = name.strip(";")
            stack[stack_index].append(name)
            if "protein" in name:
                num = int(name.split("_")[1])
                seq = protein_seqs[num]
                tree[name] = seq
            else:
                seq = node_dict[name]
                tree[name] = Node(seq)
            i = newick.find(":", i) + 1

        else:
            i = i + 1
        print(stack)
    return list(tree.values())

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
    output = build_tree(seqs, test_hom)
    for i in output:
        if isinstance(i, str):
            print(i)
        else:
            print(f"sequence:{i.sequence}, children:{i.children}, memory:{i}")
