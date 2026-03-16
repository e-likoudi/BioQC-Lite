# BioQC-Lite

BioQC-Lite is a lightweight, unified Quality Control tool designed for rapid analysis of biological sequences. Whether you are working with DNA, RNA, or Proteins, this tool provides essential assembly metrics and compositional insights in seconds.



## Features

- **Multi-Omics Support**: Automatically detects sequence types (DNA, RNA, or Protein)
- **Unified Processing**: Handles both FASTA and FASTQ formats
- **Assembly Metrics**: Calculates N50, Max/Min lengths, and total sequence counts
- **Compositional Analysis**: Provides GC Content distribution for nucleotides
- **Complexity Profiling**: Generates the top 10 most frequent k-mers
- **Visual Reports**: Automatically generates distribution plots using `matplotlib`



## Installation

You can install BioQC-Lite directly from PyPI:

```bash
  pip install bioqc-lite
```
    
## Usage

1. To initialise, open your terminal in a new folder and run `bioqc`.
    
    This will automatically create a `seq_files/` directory.

2. Add your data by placing the sequence files (.fasta, .fastq, .fq) inside the `seq_files/` folder.

3. To analyze your data, run the `bioqc` command again.

4. Open `reports/summary.txt` for detailed statistics.

    View `reports/sequence_analysis_report.png` for visual distributions.
    



## Screenshots

- Example of `summary.txt`

<img width="263" height="526" alt="Στιγμιότυπο οθόνης 2026-03-17, 1 23 29 πμ" src="https://github.com/user-attachments/assets/dca56189-9273-4e95-951a-fc6c2785deb4" />



- Example of distribution plots

<img width="640" height="480" alt="sequence_analysis_report" src="https://github.com/user-attachments/assets/ba38edb4-0adb-426a-be50-f38819775d4d" />


## Scientific Metrics

- **N50 Score**

The N50 is defined as the sequence length $L_i$ such that:

$$\sum_{j=1}^{i} L_j \ge 0.5 \times \sum_{j=1}^{n} L_j$$

It is a key metric for assessing the quality of genomic assemblies.


- **GC Content**

For DNA and RNA, the tool calculates the fraction of Guanine (G) and Cytosine (C) bases, which is critical for understanding genomic stability and gene density.
## Project Structure

```bash
BioQC-Lite/
├── src/bioqc/
│   ├── __init__.py
│   └── main.py
├── reports/          # Auto-generated reports
└── seq_files/        # User-provided sequence data
```
## 💭 Feedback and Contributing

Feel free to fork, report issues, or submit pull requests! You can also use the discussion tab for any related questions:
- https://github.com/e-likoudi/BioQC-Lite/discussions
