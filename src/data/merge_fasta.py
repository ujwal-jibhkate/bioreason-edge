import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm  # Added for progress visualization
import logging

# --- Configuration ---
RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
FILES = {
    "historical": RAW_DIR / "lanl_2010_2019.fasta",
    "recent": RAW_DIR / "lanl_2020_2025.fasta"
}
OUTPUTS = {
    "master": PROCESSED_DIR / "master_hiv_pol.fasta",
    "hyena": PROCESSED_DIR / "hyena_pretrain.txt"
}

# Thresholds
MIN_BP = 800
MAX_AMBIGUITY = 0.10  # 10%

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_accession(header_str: str) -> str:
    """
    Extracts Accession from header: >Subtype.Country.Year.Name.Accession
    Example: >B.FR.1983.HXB2.K03455 -> K03455
    """
    try:
        # Split by dot, take the last element
        parts = header_str.split('.')
        return parts[-1].strip()
    except IndexError:
        logger.warning(f"Malformed header found: {header_str}")
        return ""

def validate_and_clean(record: SeqRecord) -> SeqRecord:
    """
    Cleans sequence (removes gaps) and validates length/ambiguity.
    Returns None if validation fails.
    """
    # 1. Clean: Convert to upper case and remove alignment gaps ('-')
    # Critical for removing the '-------' seen in aligned inputs
    raw_seq_str = str(record.seq).upper().replace("-", "")
    
    # 2. Length Check
    if len(raw_seq_str) < MIN_BP:
        return None

    # 3. Ambiguity Check
    # Calculate non-standard bases (anything not A, C, G, T)
    standard_bases = set("ACGT")
    ambiguous_count = sum(1 for base in raw_seq_str if base not in standard_bases)
    ambiguity_ratio = ambiguous_count / len(raw_seq_str)

    if ambiguity_ratio > MAX_AMBIGUITY:
        return None

    # Return new record with cleaned sequence
    return SeqRecord(
        Seq(raw_seq_str),
        id=record.id,
        description=record.description,
        name=record.name
    )

def process_file(filepath, unique_sequences, stats_key_total, stats_key_duplicates=None):
    """
    Helper function to process a single FASTA file.
    """
    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        sys.exit(1)

    logger.info(f"Processing {filepath.name}...")
    
    # We use tqdm to show a progress bar. 
    # Since we don't know total lines beforehand, it will show iterations/sec.
    count = 0
    overwritten = 0
    failed = 0
    
    with open(filepath, "r") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc=filepath.name, unit="seq"):
            count += 1
            clean_record = validate_and_clean(record)
            
            if clean_record:
                acc = get_accession(record.description)
                if acc:
                    if acc in unique_sequences and stats_key_duplicates:
                        overwritten += 1
                    unique_sequences[acc] = clean_record
            else:
                failed += 1
                
    return count, overwritten, failed

def main():
    # Ensure output directory exists
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Dictionary to store unique sequences: {Accession: SeqRecord}
    unique_sequences = {}
    
    stats = {
        "batch1_total": 0, "batch1_failed": 0,
        "batch2_total": 0, "batch2_failed": 0,
        "duplicates_overwritten": 0
    }

    # --- Pass 1: Historical Data (2010-2019) ---
    b1_total, _, b1_failed = process_file(
        FILES["historical"], unique_sequences, "batch1_total"
    )
    stats["batch1_total"] = b1_total
    stats["batch1_failed"] = b1_failed

    # --- Pass 2: Recent Data (2020-2025) ---
    # Recent overrides Old
    b2_total, b2_overwritten, b2_failed = process_file(
        FILES["recent"], unique_sequences, "batch2_total", stats_key_duplicates=True
    )
    stats["batch2_total"] = b2_total
    stats["batch2_failed"] = b2_failed
    stats["duplicates_overwritten"] = b2_overwritten

    # --- Output Generation ---
    final_records = list(unique_sequences.values())
    
    # 1. Write Master FASTA
    logger.info(f"Writing {len(final_records)} sequences to {OUTPUTS['master']}")
    SeqIO.write(final_records, OUTPUTS["master"], "fasta")

    # 2. Write Hyena Pretrain (Raw Text)
    logger.info(f"Writing raw sequences to {OUTPUTS['hyena']}")
    with open(OUTPUTS["hyena"], "w") as f:
        for record in final_records:
            f.write(str(record.seq) + "\n")

    # --- Final Summary ---
    print("\n" + "="*40)
    print(" PIPELINE SUMMARY")
    print("="*40)
    print(f"Batch 1 (Historical): {stats['batch1_total']} seqs (Failed val: {stats['batch1_failed']})")
    print(f"Batch 2 (Recent):     {stats['batch2_total']} seqs (Failed val: {stats['batch2_failed']})")
    print(f"Duplicates Resolved:  {stats['duplicates_overwritten']} (Recent ver. kept)")
    print("-" * 40)
    print(f"FINAL TOTAL SEQUENCES:      {len(final_records)}")
    print(f"HXB2 Status (K03455):       {'Found' if 'K03455' in unique_sequences else 'Missing'}")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()