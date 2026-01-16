import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import logging

# --- Configuration ---
RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
FILES = {
    "historical": RAW_DIR / "lanl_2010_2019.fasta",
    "recent": RAW_DIR / "lanl_2020_2025.fasta"
}
OUTPUTS = {
    "master": PROCESSED_DIR / "master_hiv_genome.fasta", # Renamed for clarity
    "hyena": PROCESSED_DIR / "hyena_pretrain_genome.txt"
}

# Thresholds for WHOLE GENOME
# HIV-1 Genome is approx 9.7kb. We set 8000 to allow for small terminal deletions.
MIN_BP = 8000 
MAX_AMBIGUITY = 0.10  # 10%

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_accession(header_str: str) -> str:
    """
    Extracts Accession from header: >Subtype.Country.Year.Name.Accession
    """
    try:
        parts = header_str.split('.')
        return parts[-1].strip()
    except IndexError:
        logger.warning(f"Malformed header found: {header_str}")
        return ""

def validate_and_clean(record: SeqRecord) -> SeqRecord:
    """
    Cleans sequence (removes gaps) and validates length/ambiguity.
    """
    # 1. Clean: Upper case and remove alignment gaps
    raw_seq_str = str(record.seq).upper().replace("-", "")
    
    # 2. Length Check (Updated for Genome)
    if len(raw_seq_str) < MIN_BP:
        return None

    # 3. Ambiguity Check
    standard_bases = set("ACGT")
    ambiguous_count = sum(1 for base in raw_seq_str if base not in standard_bases)
    ambiguity_ratio = ambiguous_count / len(raw_seq_str)

    if ambiguity_ratio > MAX_AMBIGUITY:
        return None

    return SeqRecord(
        Seq(raw_seq_str),
        id=record.id,
        description=record.description,
        name=record.name
    )

def process_file(filepath, unique_sequences, stats_key_total, stats_key_duplicates=None):
    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        sys.exit(1)

    logger.info(f"Processing {filepath.name}...")
    
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
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    unique_sequences = {}
    
    stats = {
        "batch1_total": 0, "batch1_failed": 0,
        "batch2_total": 0, "batch2_failed": 0,
        "duplicates_overwritten": 0
    }

    # --- Pass 1: Historical ---
    b1_total, _, b1_failed = process_file(
        FILES["historical"], unique_sequences, "batch1_total"
    )
    stats["batch1_total"] = b1_total
    stats["batch1_failed"] = b1_failed

    # --- Pass 2: Recent ---
    b2_total, b2_overwritten, b2_failed = process_file(
        FILES["recent"], unique_sequences, "batch2_total", stats_key_duplicates=True
    )
    stats["batch2_total"] = b2_total
    stats["batch2_failed"] = b2_failed
    stats["duplicates_overwritten"] = b2_overwritten

    # --- Output ---
    final_records = list(unique_sequences.values())
    
    logger.info(f"Writing {len(final_records)} sequences to {OUTPUTS['master']}")
    SeqIO.write(final_records, OUTPUTS["master"], "fasta")

    logger.info(f"Writing raw sequences to {OUTPUTS['hyena']}")
    with open(OUTPUTS["hyena"], "w") as f:
        for record in final_records:
            f.write(str(record.seq) + "\n")

    # --- Summary ---
    print("\n" + "="*40)
    print(" PIPELINE SUMMARY (WHOLE GENOME)")
    print("="*40)
    print(f"Batch 1 Input: {stats['batch1_total']} (Dropped <8000bp/Qual: {stats['batch1_failed']})")
    print(f"Batch 2 Input: {stats['batch2_total']} (Dropped <8000bp/Qual: {stats['batch2_failed']})")
    print(f"Duplicates Resolved: {stats['duplicates_overwritten']}")
    print("-" * 40)
    print(f"FINAL TOTAL SEQUENCES: {len(final_records)}")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()