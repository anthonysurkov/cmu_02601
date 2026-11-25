import pathlib

import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer


MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
REPRESENTATION_LAYER = 33
BATCH_SIZE = 2
MAX_LENGTH = 1022


def get_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_model(device: torch.device):
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)

    model.eval()
    model.to(device)

    return tokenizer, model


def extract_embeddings(
    sequences: list[str],
    tokenizer,
    model,
    device: torch.device,
) -> list[list[float]]:
    embeddings = []

    for start in range(0, len(sequences), BATCH_SIZE):
        batch = sequences[start:start + BATCH_SIZE]

        tokens = tokenizer(
            batch,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=MAX_LENGTH + 2,
        ).to(device)

        with torch.inference_mode():
            outputs = model(
                **tokens,
                output_hidden_states=True,
            )

        hidden = outputs.hidden_states[REPRESENTATION_LAYER]

        attention_mask = tokens["attention_mask"]

        for i in range(len(batch)):
            sequence_length = attention_mask[i].sum().item()

            residue_embeddings = hidden[i, 1:sequence_length - 1]

            mean_embedding = residue_embeddings.mean(dim=0)

            embeddings.append(
                mean_embedding.cpu().tolist()
            )

        print(f"Processed {min(start + BATCH_SIZE, len(sequences))}"
              f"/{len(sequences)} sequences")

    return embeddings


def process_split(
    csv_path: pathlib.Path,
    output_path: pathlib.Path,
    tokenizer,
    model,
    device: torch.device,
) -> None:
    df = pd.read_csv(csv_path)

    sequences = df["sequence"].tolist()

    embeddings = extract_embeddings(
        sequences,
        tokenizer,
        model,
        device,
    )

    embedding_df = pd.DataFrame(
        embeddings,
        columns=[
            f"embedding_{i}" for i in range(len(embeddings[0]))
        ],
    )

    result = pd.concat(
        [
            df[["ProteinID", "Tm"]].reset_index(drop=True),
            embedding_df,
        ],
        axis=1,
    )

    result.to_csv(output_path, index=False)

    print(f"Saved embeddings to {output_path}")


def main():
    data_dir = pathlib.Path("data/tmprot")
    output_dir = pathlib.Path("data/tmprot/embeddings")

    output_dir.mkdir(parents=True, exist_ok=True)

    device = get_device()
    print(f"Using device: {device}")

    tokenizer, model = load_model(device)

    for split in ["train", "val", "test"]:
        process_split(
            csv_path=data_dir / f"{split}.csv",
            output_path=output_dir / f"{split}_embeddings.csv",
            tokenizer=tokenizer,
            model=model,
            device=device,
        )


if __name__ == "__main__":
    main()

