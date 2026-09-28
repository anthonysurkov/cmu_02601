# 02601 CMU: Semester Project
### Proposal:
https://docs.google.com/document/d/1CBZraYHg-huKTE3e-L1dlmMjN-Q35B27imOgkhhG3mA/edit?usp=sharing
Anthony Surkov, Anushka Shome, Gauri Chasia, Vaishnavi Venuturimilli 

### Project management:
(0) Basic infrastructure (Not assigned)
- References (e.g. canonical protein reference), caching, path management
(0a) Visualization (Not assigned)
- Protein rotation gifs with proposed mutations highlighted (most protein sim software can do this)

(1) Deep mutational scanning via TemStaPro model (Not assigned)
- Inputs: str of canonical BglB enzyme sequence
- Outputs: dict[int, str] of all mutation positions and new amino acids that are predicted to increase thermostability
(2) BLAST API wrapper (Not assigned)
- Inputs: str of canonical BglB enzyme sequence
- Outputs: list[str] of homologous BglB enzyme sequences
(3a) IQ-TREE ModelFinder wrapper (Not assigned)
- Inputs: list[str] of homologous BglB enzyme sequences
- Outputs: (type unknown) amino-acid substitution model for use with ASR
(3b) ASR (Not assigned)
- Inputs: list[str] of homologous BglB enzyme sequences, amino-acid substitution model
- Outputs: list[Node | str] of extant and ancestral BglB enzyme sequences (flattened from tree)
(4) Alignment (Anthony - not started)
- Takes: list[str] of extant and ancestral BglB enzyme sequences
- Outputs: list[str] of aligned sequences, list[float] of agreement between residues across aligned sequences
(4a) Alignment decision algorithm (Not assigned)
- Takes: list[str] of aligned sequences, list[float] of agreement between residues across aligned sequences
- Outputs: dict[int, str] of all positions strongly conserved in ASR (residue agreement within some parameter theta)
(5) Final decision algorithm (Not assigned)
- Takes: dict[int, str] of all mutation positions and amino acids, dict[int, str] of conserved positions
- Outputs: dict[int, str] of cross-referenced mutations thought to be viable for both thermostability and ASR

(6) Feature modules (Not assigned)
- Takes: Downloaded CSV from D2D database OR proposed mutations from (5)
- Outputs: Feature cache of ProtBERT descriptors per protein sequence provided
(7) D2D BglB database truth-modules (Not assigned)
- Takes: Feature cache from (6)
- Outputs: normalized torch.utils.data.Dataset objects for model training, canonical train/test splits
(8) Training & eval modules (Not assigned)
- Takes: torch.utils.data.Dataset objects
- Outputs: XGBoost (or other) models (.pt)
(9) Plotting modules (Not assigned)
- Takes: y_val, y vectors from Downloaded CSV, (8)'s model predictions
- Outputs: Scatterplots of model performance; model diagnostics
