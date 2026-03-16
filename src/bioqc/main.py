from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from matplotlib import pyplot as plt
from collections import Counter

SEQ_PATH = Path("seq_files")
REPORT_PATH = Path("reports")

# Load sequence files from the specified directory and categorize them by format
def load_seq_files() -> list:
    fasta_files = []
    fastq_files = []
    if not SEQ_PATH.exists():
        return [], []

    for filename in SEQ_PATH.iterdir():
        if filename.suffix == ".fasta":
            fasta_files.append(filename)
        elif filename.suffix == ".fastq" or filename.suffix == ".fq":
            fastq_files.append(filename)
    return fasta_files, fastq_files

# Determine if a sequence is DNA, RNA, or Protein based on its characters
def detect_sequence_type(sequence: str) -> str:
    capsequence = sequence.upper()
    seq_set = set(capsequence) 
    nucleotides = {'A', 'C', 'G', 'T', 'U', 'N'}
    
    if not seq_set.issubset(nucleotides):
        return "Protein"
    if 'U' in seq_set:
        return "RNA"
    return "DNA"

# Generate k-mers from a given sequence. If the sequence is shorter than k, return an empty list.
def generate_kmers(sequence: str, k: int):
    sec_lenght = len(sequence)
    i = 0
    if sec_lenght < k:
        print(f"Sequence length {sec_lenght} is shorter than k={k}. Skipping k-mer generation.")
        return []
    
    for i in range(sec_lenght - k + 1):
        yield sequence[i:i+k]

# Process each sequence file, updating statistics and k-mer counts. This function handles both FASTA and FASTQ formats.
def process_sequence_files(file_list: list, format_type: str, stats: dict, kmer_stats: Counter, k: int) -> None:
    for filename in file_list:
        print(f"Processing {format_type.upper()} file: {filename}")
        filepath = SEQ_PATH / filename
        for record in SeqIO.parse(filepath, format_type):
            seq_str = str(record.seq)
            seq_type = detect_sequence_type(seq_str)
            stats[seq_type]["lengths"].append(len(seq_str))

            if seq_type in ["DNA", "RNA"]:
                stats[seq_type]["gc"].append(gc_fraction(seq_str))

            kmers = generate_kmers(seq_str, k)
            kmer_stats.update(kmers)

# Calculate the N50 value from a list of sequence lengths.
def calculate_n50(lengths: list):
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    return 0

# Plot the distribution of sequence lengths and GC content using Matplotlib.    
def plot_length_distribution(stats: dict):
    all_lengths = stats["DNA"]["lengths"] + stats["RNA"]["lengths"] + stats["Protein"]["lengths"]
    nucleotides_gc = stats["DNA"]["gc"] + stats["RNA"]["gc"]

    plt.subplot(1, 2, 1)
    plt.hist(all_lengths, bins=20, color='skyblue', edgecolor='black')
    plt.title("Sequence Length Distribution")
    plt.xlabel("Length (bp/aa)")
    plt.ylabel("Frequency")
    plt.grid(axis='y', alpha=0.5)

    if nucleotides_gc:
        plt.subplot(1, 2, 2)
        plt.hist(nucleotides_gc, bins=20, range=(0, 1), color='lightgreen', edgecolor='black', alpha=0.7)
        plt.title('GC Content Distribution\n (Nucleotides Only)')
        plt.xlabel('GC %')
        plt.ylabel('Frequency')
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0%', '20%', '40%', '60%', '80%', '100%'])
        plt.grid(axis='y', alpha=0.5)
    else:
        plt.subplot(1, 2, 2)
        plt.text(0.5, 0.5, 'No Nucleotide sequences found', horizontalalignment='center', verticalalignment='center')
        plt.title('GC Content Distribution\n (No Data)')
        plt.axis('off')

    plt.tight_layout()
    output_filename = "sequence_analysis_report.png"
    plt.savefig(REPORT_PATH / output_filename) 
    plt.close()

# Save a summary of the analysis, including total sequences, N50 values, and top k-mers, to a text file.
def save_summary(stats: dict, kmer_stats: Counter):
    all_lengths = stats["DNA"]["lengths"] + stats["RNA"]["lengths"] + stats["Protein"]["lengths"]
    
    with open(REPORT_PATH / "summary.txt", "w") as f:
        f.write("Sequence Type Summary:\n\n")
        f.write(f"Total Sequences Analyzed: {len(all_lengths)}\n")
        f.write(f"Unique K-mers found:      {len(kmer_stats)}\n")
        
        for seq_type in ["DNA", "RNA", "Protein"]:
            lengths = stats[seq_type]["lengths"]
            gc_values = stats[seq_type]["gc"]
            
            if lengths:
                f.write(f"\nType: {seq_type}\n")
                f.write(f"  - Count:  {len(lengths)}\n")
                f.write(f"  - N50:    {calculate_n50(lengths)} bp/aa\n")
                f.write(f"  - Max:    {max(lengths)} | Min: {min(lengths)}\n\n")
                
                if gc_values:
                    avg_gc = sum(gc_values) / len(gc_values)
                    f.write(f"  - Avg GC: {avg_gc:.2%}\n")

        f.write(f"Global N50: {calculate_n50(all_lengths)}\n")
        if kmer_stats:
            f.write(f"Top 10 most frequent k-mers:\n")
            for kmer, count in kmer_stats.most_common(10):
                f.write(f"  - {kmer}: {count}\n")
        else:
            f.write("No k-mers generated.\n")

def main():
    SEQ_PATH.mkdir(parents=True, exist_ok=True)
    fasta_files, fastq_files = load_seq_files()
    print("FASTA files found:", len(fasta_files))
    print("FASTQ files found:", len(fastq_files))
    if not fasta_files and not fastq_files:
        print(f"No FASTA or FASTQ files found in the directory. Please place your files in the {SEQ_PATH} folder and run again")
    REPORT_PATH.mkdir(parents=True, exist_ok=True)

    k = input("Enter the value of k for k-mer generation: ")
    try:
        k = int(k)
        if k <= 0:
            print("k must be a positive integer. Exiting.")
            return
    except ValueError:
        print("Invalid input for k. Please enter a positive integer. Exiting.")
        return

    kmer_stats = Counter()
    stats = {
        "DNA": {"lengths": [], "gc": []},
        "RNA": {"lengths": [], "gc": []},
        "Protein": {"lengths": [], "gc": []}
    }

    process_sequence_files(fasta_files, "fasta", stats, kmer_stats, k)
    print(f"Processed {len(fasta_files)} FASTA files.")
    process_sequence_files(fastq_files, "fastq", stats, kmer_stats, k)
    print(f"Processed {len(fastq_files)} FASTQ files.")

    plot_length_distribution(stats)
    print(f"Length distribution plot saved to {REPORT_PATH}sequence_analysis_report.png")

    save_summary(stats, kmer_stats)
    print(f"Summary saved to {REPORT_PATH}summary.txt")
    
if __name__ == "__main__":
    main()