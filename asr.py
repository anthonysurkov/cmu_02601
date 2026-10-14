'''
ASR (part 3b)
Inputs (assumed, will change if needed):
  - list[str] containing only the sequences, no labels
  - str with model type (eg. "LG")

Outputs (will change depending on next step):
  - Returns a list of mixed types: type Node for internal nodes + root, and
      type str for extant proteins
  - Folder with iqtree output, labeled "iqtree_output"
  - FASTA file with extant proteins labeled "homologous_bglb_sequences.fasta"

Requires iqtree2 to be installed in environment:
    conda install -c bioconda iqtree=2
'''


from pathlib import Path
import subprocess
from node import Node


fasta_path = Path("homologous_bglb_sequences.fasta")
output_dir = Path("iqtree_output")


def create_fasta(sequences: list[str]) -> None:

    '''
    Takes in the list of extant protein sequences from part 3a
    Writes them into a FASTA file with arbitrary labels (protein_[num])
    Used for running IQTree in the next step
    '''

    with fasta_path.open("w") as f:
        for i in range(len(sequences)):
            f.write(f">protein_{i}\n")
            f.write(f"{sequences[i]}\n")


def run_iqtree(model: str) -> None:

    '''
    Takes in the model type from the part 3a
    Output directory iqtree_output is created and iqtree output is stored in
        there
    Note: Will fail if iqtree_output directory already exists, so make sure to
          rm -r if testing multiple times
    '''

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / "bglb_asr"

    command = ["iqtree2", "-s", fasta_path, "-st", "AA", "-m", model, "-asr",
        "-pre", str(prefix)
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


def extract_sequences() -> dict[str, str]:

    '''
    Extract sequences for each node of tree from the iqtree state file
    Retuns a dictionary, where the key is the node name and the value is the
        reconstructed sequence at that node
    '''

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


def build_tree(node_dict: dict[str, str],
    protein_seqs: list[str]
) -> list[Node | str]:

    '''
    Takes in the node dictionary made by extract_sequences
    Also takes in the extant sequences given as input to ASR
    Returns a mix of Node and strings, representing ancestral nodes made by
      IQTree and extant sequences, respectively

    Deciphers the Newick Format from the tree file in iqtree_output/:
       - Data structures:
           - tree: dictionary keeping track of the name of each node (Node[num]
             or protein_num) and its corresponding item (Node or str)
           - stack: dictionary keeping track of the sibling groups, where the
             key is a number corresponding to each sibling group and the values
             are lists of all the nodes in each group.
               - sibling group key index is kept track of through stack_index
       - If there is an opening parenthesiss "(", it means a new sibling group
         is starting. Level of stack is incremented, and an item is appended to
         the tree dictionary.
       - If there is a closing parenthesis ")", it means the current sibling
         group has ended, and the protein/Node listed immediately after is the
         direct parent of all the siblings in the current group. A new Node
         object is created for the parent with all the siblings listed as
         children, and is added to the tree dictionary  Level of stack is
         decremented, and the parent node is added to that sibling group in the
         stack.
       - If there is a comma "," and it is followed by alphabet characters, it
         must correspond to another protein/Node that is in the current sibling
         group. An item is appended to the tree dictionary and added to the
         current sibling group in the stack.
    '''

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
            else: # Done (will always end with the parent node name next to ')')
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

    return list(tree.values())


if __name__ == "__main__":
    # Test placeholders for now, will get rid of later
    test_hom = ["MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNTD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNTE",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVSTD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVSTE",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNAD",
        "MKKLAVVLAGASVAGSAFAAPVNTTQKLEQALKKAGVEVVNAE"
    ]
    model_type = "LG"

    create_fasta(test_hom)
    run_iqtree(model_type)
    seqs = extract_sequences()
    output = build_tree(seqs, test_hom)
