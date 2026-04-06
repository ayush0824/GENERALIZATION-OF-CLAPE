# CLAPE-RNA: Predicting RNA-Protein Binding Sites

**Author:** Ayush Chhangani  
**Thesis:** *Predict RNA-Protein Binding Sites Using the CLAPE Framework* 
**Institution:** The Pennsylvania State University, J. Jeffrey and Ann Marie Fox Graduate School

##  Key Features
- **Transformer-based Embeddings:** Utilizes ProtBert for high-dimensional protein representation.
- **Contrastive Learning:** Employs a locality-aware contrastive task to identify binding residues.
- **Adaptive Positional Embeddings:** Specifically tuned for RNA sequence context dependence.

## Data & Preprocessing

### Custom 3-Line Format
Our pipeline uses a custom 3-line format to accommodate experimental masks:
1.  **Header** (e.g., `>7zod_A`)
2.  **Sequence** (e.g., `AVQ...`)
3.  **Mask/Label** (Binary or Probability)

**Important:** To convert our custom format to standard FASTA for external tools, use the provided utility:
```bash
python GetFilteredProteinsFromTestSet.py <input_3line.txt> <output.fasta>
```

Rigorous Scientific Integrity (PDB-30)
In bioinformatics, your model is only as good as your validation. Unlike many tools that allow data leakage, CLAPE-RNA is built on a High-Fidelity Evaluation Protocol:

Strict 30% Identity Filter: We use psi-cd-hit to ensure the test set contains no sequences with >30% homology to the training data.

40% Coverage Guardrail: Validating that global protein architecture doesn't bias the local site prediction.

Custom 3-Line Processing: Native support for experimental masks and probability labels, ensuring data integrity from raw file to final output.

To build the filtered-protein set

To run CD-HIT pipeline:

	1) sed 's/^>/>TRAIN_/' train.fasta > train_tagged.fasta

	2) sed 's/^>/>TEST_/' test.fasta > test_tagged.fasta
	
	3) cat train_tagged.fasta test_tagged.fasta > combined.fasta

	4) psi-cd-hit.pl -i combined.fasta -o combined_pdb30.fasta -c 0.3

	5) 	python3 parse-pool-cluster-cd-hit-proteins.py \
 		 --clstr combined_pdb30.fasta.clstr \ 
  		--test RNA_1500_Test.txt \
 		 --out T-2025_PDB30_Final.txt

Here, in last command, all files are taken as an input dynamically. 


To build files for clape predict:

 python3 GetFilteredProteinsFromTestSet.py  T-2025_PDB30_Final.txt T-PDB30-final-filtered.fasta

Performance & Results
Validated on BioLip and ENCODE datasets: Curated using large-scale BioLip, eCLIP-seq and PDB datasets.

State-of-the-Art Metrics: 
Optimized for AUPRC, AUROC, and F1 scores to provide researchers with reliable, actionable insights.

Citation
If CLAPE-RNA accelerates your research, please cite our work:

The Thesis:

Chhangani, A. (2025). Predict RNA-Protein Binding Sites Using the CLAPE Framework. Master's Thesis, The J. Jeffrey and Ann Marie Fox Graduate School, The Pennsylvania State University.

The Foundation:

Liu, Y., & Tian, B. (2024). Protein–DNA binding sites prediction based on pre-trained protein language model and contrastive learning. Briefings in Bioinformatics.

Maintained by Ayush Chhangani Computer Science & Engineering, Penn State University
