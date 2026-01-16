import sys
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import logging

# --- Configuration ---
# We assume the balanced fasta is in the processed directory
PROCESSED_DIR = Path("data/processed")

INPUT_FASTA = PROCESSED_DIR / "balanced_hiv_genome.fasta"
OUTPUT_TXT = PROCESSED_DIR / "hyena_pretrain_genome_balanced.txt" 

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    if not INPUT_FASTA.exists():
        logger.error(f"Input file not found: {INPUT_FASTA}")
        print("Please ensure you ran the CD-HIT command successfully in the 'data/processed' directory.")
        sys.exit(1)

    logger.info(f"Reading from: {INPUT_FASTA}")
    logger.info(f"Writing to:   {OUTPUT_TXT}")

    count = 0
    
    # We open the output file in write mode
    with open(OUTPUT_TXT, "w") as f_out:
        # SeqIO handles multi-line FASTA records automatically
        for record in tqdm(SeqIO.parse(INPUT_FASTA, "fasta"), desc="Extracting sequences", unit="seq"):
            
            # 1. Get the sequence
            # CD-HIT might wrap lines, but SeqIO.parse joins them back into one string.
            # We explicitly uppercase and remove gaps just to be 100% safe, 
            # though your previous script likely already did this.
            seq_str = str(record.seq).upper().replace("-", "")
            
            # 2. Write to text file (One sequence per line)
            f_out.write(seq_str + "\n")
            
            count += 1

    # --- Summary ---
    print("\n" + "="*40)
    print(" TXT GENERATION COMPLETE")
    print("="*40)
    print(f"Source FASTA: {INPUT_FASTA.name}")
    print(f"Output TXT:   {OUTPUT_TXT.name}")
    print("-" * 40)
    print(f"TOTAL SEQUENCES WRITTEN: {count}")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()