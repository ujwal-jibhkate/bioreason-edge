import os
import json
import time
import requests
import logging
from typing import List, Dict, Any, Tuple
from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

# --- Configuration ---
CONFIG = {
    "INPUT_FILE": "data/processed/balanced_hiv_genome.fasta",
    "OUTPUT_FILE": "data/processed/train_instruction.jsonl",
    "API_URL": "https://hivdb.stanford.edu/graphql",
    "BATCH_SIZE": 5, 
    "RATE_LIMIT_DELAY": 1.0,  
    "MAX_RETRIES": 3
}

# --- Logging Setup ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("dataset_generation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# --- CORRECTED GRAPHQL QUERY ---
# Sourced directly from your working example.
# We explicitly request genes [CA, PR, RT, IN] to ensure consistent output.
SIERRA_QUERY = """
query ($sequences: [UnalignedSequenceInput]!) {
  sequenceAnalysis(sequences: $sequences) {
    
    inputSequence { header }

    # Mutations (Source of Truth)
    alignedGeneSequences(includeGenes: [CA, PR, RT, IN]) {
      gene { name }
      mutations {
        text
        shortText
        primaryType
      }
    }

    # Drug Resistance Scores (Returns a LIST of gene reports)
    drugResistance(includeGenes: [CA, PR, RT, IN]) {
      gene { name }
      drugScores {
        drugClass { name }
        drug { name displayAbbr }
        score
        text
        level
      }
    }
  }
}
"""

def get_session():
    """Creates a requests session with robust retry logic."""
    session = requests.Session()
    retry = Retry(
        total=CONFIG["MAX_RETRIES"],
        read=CONFIG["MAX_RETRIES"],
        connect=CONFIG["MAX_RETRIES"],
        backoff_factor=2, 
        status_forcelist=[500, 502, 503, 504],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session

def count_processed_sequences(filepath: str) -> int:
    """Counts lines in JSONL to support resume capability."""
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r', encoding='utf-8') as f:
        return sum(1 for _ in f)

def read_fasta_in_batches(filepath: str, batch_size: int, skip_count: int):
    """Yields batches of (header, sequence) tuples."""
    batch = []
    current_header = None
    current_seq = []
    processed_count = 0
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                if current_header:
                    processed_count += 1
                    if processed_count > skip_count:
                        batch.append((current_header, "".join(current_seq)))
                        if len(batch) == batch_size:
                            yield batch
                            batch = []
                
                current_header = line[1:] 
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            processed_count += 1
            if processed_count > skip_count:
                batch.append((current_header, "".join(current_seq)))
                if batch:
                    yield batch

def clean_sequence(seq: str) -> str:
    """
    Sanitizes DNA string for the API.
    API fails if there are gaps/newlines/spaces.
    """
    return seq.upper().replace("-", "").replace(" ", "").replace("\n", "").strip()

def format_clinical_report(analysis: Dict[str, Any]) -> str:
    """Formats the JSON analysis into a structured report."""
    if not analysis:
        return "Clinical Report: Analysis failed."
    
    # 1. Process Mutations (From alignedGeneSequences LIST)
    mutations_text = []
    gene_seqs = analysis.get('alignedGeneSequences', [])
    
    if not gene_seqs:
         return "Clinical Report: No HIV-1 Polymerase genes (PR, RT, IN) were detected in this sequence."

    for gene_obj in gene_seqs:
        gene_name = gene_obj['gene']['name']
        muts = [m['text'] for m in gene_obj.get('mutations', [])]
        if muts:
            mutations_text.append(f"{gene_name}: {', '.join(muts)}")
    
    mutations_section = "\n".join(mutations_text) if mutations_text else "No major surveillance mutations detected."

    # 2. Process Drug Resistance (From drugResistance LIST)
    dr_list = analysis.get('drugResistance', [])
    if dr_list is None: dr_list = []
    
    drug_classes = {}
    
    # Iterate through the LIST of gene reports
    for gene_report in dr_list:
        for score in gene_report.get('drugScores', []):
            d_class = score['drugClass']['name']
            d_name = score['drug']['name']
            d_abbr = score['drug']['displayAbbr']
            d_text = score['text'] # e.g. "High-Level Resistance"
            
            if d_class not in drug_classes:
                drug_classes[d_class] = []
            
            # Create entry string
            entry = f"- {d_name} ({d_abbr}): {d_text}"
            
            # Deduplicate (sometimes multiple genes report same drug)
            if entry not in drug_classes[d_class]:
                drug_classes[d_class].append(entry)

    profile_text = []
    for d_class, drugs in drug_classes.items():
        profile_text.append(f"**{d_class}**")
        profile_text.extend(drugs)
        profile_text.append("") 

    profile_section = "\n".join(profile_text).strip()

    return f"""### Clinical HIV Drug Resistance Report

**Genetic Mutations Identified:**
{mutations_section}

**Resistance Profile:**
{profile_section}"""

def generate_reasoning_trace(analysis: Dict[str, Any]) -> str:
    """Generates a synthetic 'Chain-of-Thought' reasoning trace."""
    steps = []
    steps.append("Step 1: Analyzing the viral genome sequence for known resistance motifs (PR, RT, IN).")

    # Collect mutations
    found_muts = []
    gene_seqs = analysis.get('alignedGeneSequences', [])
    if gene_seqs:
        for gene_obj in gene_seqs:
            g_name = gene_obj['gene']['name']
            for m in gene_obj.get('mutations', []):
                 found_muts.append(f"{g_name}:{m['text']}")
    
    if found_muts:
        display_muts = found_muts[:8]
        steps.append(f"Step 2: Identified surveillance mutations: {', '.join(display_muts)}" + ("..." if len(found_muts)>8 else "."))
    else:
        steps.append("Step 2: No major drug resistance mutations identified.")

    steps.append("Step 3: Cross-referencing identified mutations with the Stanford HIVDB scoring algorithm.")

    # Check Resistance (Iterating the LIST)
    high_res_drugs = []
    dr_list = analysis.get('drugResistance', [])
    if dr_list:
        for gene_report in dr_list:
            for score in gene_report.get('drugScores', []):
                # Check for High or Intermediate
                if "High" in score['text'] or "Intermediate" in score['text']:
                    high_res_drugs.append(score['drug']['displayAbbr'])
    
    # Deduplicate
    high_res_drugs = list(set(high_res_drugs))

    if high_res_drugs:
        steps.append(f"Step 4: Significant resistance penalties detected for: {', '.join(high_res_drugs[:5])}.")
        steps.append("Step 5: Synthesizing final clinical report based on cumulative resistance scores.")
    else:
        steps.append("Step 4: No significant resistance penalties detected.")
        steps.append("Step 5: Generating standard susceptibility report.")

    return "\n".join(steps)

def process_batch(session, batch: List[Tuple[str, str]]) -> List[Dict]:
    """
    Sends a batch to Sierra API.
    Input: [ (header, sequence), ... ]
    """
    # 1. Clean Sequences & Form Payload Objects
    payload_objects = []
    for h, s in batch:
        payload_objects.append({
            "header": h,
            "sequence": clean_sequence(s)
        })
    
    # GraphQL Payload
    payload = {
        "query": SIERRA_QUERY, 
        "variables": {"sequences": payload_objects}
    }
    
    try:
        # Attempt Batch
        response = session.post(
            CONFIG["API_URL"],
            json=payload,
            timeout=60
        )
        
        # Explicit check for 400 Bad Request
        if response.status_code == 400:
            logger.error(f"API Rejected Batch. Reason: {response.text}")
            raise requests.exceptions.RequestException("Batch 400 Error")
            
        response.raise_for_status()
        
        data = response.json()
        if "errors" in data:
            logger.warning(f"GraphQL Logic Errors: {data['errors']}")
            return []

        results = data.get("data", {}).get("sequenceAnalysis", [])
        return _parse_results(batch, results)

    except (requests.exceptions.RequestException, ValueError) as e:
        logger.warning(f"Batch failed. Falling back to single processing. Error: {str(e)}")
        
        # Fallback Logic (One-by-One)
        recovered_examples = []
        for i, (header, raw_seq) in enumerate(batch):
            try:
                single_obj = payload_objects[i] # reuse cleaned object
                single_payload = {
                    "query": SIERRA_QUERY, 
                    "variables": {"sequences": [single_obj]}
                }
                
                res = session.post(CONFIG["API_URL"], json=single_payload, timeout=20)
                
                if res.status_code != 200:
                    continue
                    
                data = res.json()
                if "errors" in data or not data.get("data"):
                    continue

                analysis = data.get("data", {}).get("sequenceAnalysis", [])[0]
                parsed = _parse_results([(header, raw_seq)], [analysis])
                recovered_examples.extend(parsed)
                
            except Exception:
                continue 
        
        return recovered_examples

def _parse_results(original_batch, analysis_results):
    output_examples = []
    
    for i, analysis in enumerate(analysis_results):
        if i >= len(original_batch): break 
        
        raw_seq = original_batch[i][1]
        
        if not analysis: continue

        try:
            clinical_report = format_clinical_report(analysis)
            reasoning_trace = generate_reasoning_trace(analysis)
            
            example = {
                "messages": [
                    {
                        "role": "user",
                        "content": f"<dna_start>{raw_seq}<dna_end>\nAnalyze this HIV viral sequence and provide a clinical drug resistance report."
                    },
                    {
                        "role": "assistant",
                        "content": f"<think>\n{reasoning_trace}\n</think>\n\n{clinical_report}"
                    }
                ]
            }
            output_examples.append(example)
        except Exception as e:
            logger.error(f"Error parsing analysis for sequence {i}: {e}")
            continue

    return output_examples

def main():
    Path(os.path.dirname(CONFIG["OUTPUT_FILE"])).mkdir(parents=True, exist_ok=True)
    
    processed_count = count_processed_sequences(CONFIG["OUTPUT_FILE"])
    if processed_count > 0:
        logger.info(f"Resuming generation. Found {processed_count} existing examples.")
    else:
        logger.info("Starting new dataset generation.")

    session = get_session()
    # Estimate total based on file size or known count
    total_estimated = 3948 
    
    batch_gen = read_fasta_in_batches(CONFIG["INPUT_FILE"], CONFIG["BATCH_SIZE"], processed_count)
    pbar = tqdm(total=total_estimated, initial=processed_count, desc="Processing Genomes", unit="seq")
    
    with open(CONFIG["OUTPUT_FILE"], 'a', encoding='utf-8') as f_out:
        for batch in batch_gen:
            examples = process_batch(session, batch)
            if examples:
                for ex in examples:
                    f_out.write(json.dumps(ex) + "\n")
            
            pbar.update(len(batch))
            time.sleep(CONFIG["RATE_LIMIT_DELAY"])
            
    pbar.close()
    logger.info("Dataset generation complete.")

if __name__ == "__main__":
    main()