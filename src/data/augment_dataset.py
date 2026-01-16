import json
import random
import re
from tqdm import tqdm

# --- Configuration ---
INPUT_FILE = "data/processed/train_instruction.jsonl"
OUTPUT_FILE = "data/processed/train_instruction_augmented.jsonl"
TARGET_TOTAL_RESISTANT = 3000  # Target total count of resistant examples

# --- Genetic Mutation Knowledge Base (Multi-Subtype Support) ---
KNOWN_MUTATIONS = [
    {
        # MUTATION 1: M184V (The most common Nuke mutation)
        # Target: YMDD Motif (Tyrosine-Methionine-Aspartate-Aspartate) - Ultra Conserved
        "name": "RT: M184V",
        "gene": "Reverse Transcriptase",
        "drug_class": "NRTI",
        "drugs": ["LMV (3TC)", "FTC (FTC)"],
        "description": "High-Level Resistance",
        # We try these patterns in order. If one matches, we swap it.
        "swaps": [
            # Standard Subtype B (YMDD -> YVDD)
            ("ATGGATGAC", "GTGGATGAC"), 
            # Common Variant (YMDD -> YVDD) - Shorter anchor
            ("ATGGATGAT", "GTGGATGAT"),
        ]
    },
    {
        # MUTATION 2: K103N (The most common Non-Nuke mutation)
        "name": "RT: K103N",
        "gene": "Reverse Transcriptase",
        "drug_class": "NNRTI",
        "drugs": ["EFV (EFV)", "NVP (NVP)"],
        "description": "High-Level Resistance",
        "swaps": [
            # K103 (AAA) -> N103 (AAC) inside "Pro-Val-Lys-Lys"
            ("CCAGTAAAATTA", "CCAGTAAACTTA"), 
            ("CCAGTAAAAAAA", "CCAGTAAACAAA"),
            # Shorter anchor "Val-Lys-Lys"
            ("GTAAAATTA", "GTAAACTTA"),
        ]
    },
    {
        # MUTATION 3: D30N (Protease Inhibitor)
        "name": "PR: D30N",
        "gene": "Protease",
        "drug_class": "PI",
        "drugs": ["NFV (NFV)"],
        "description": "High-Level Resistance",
        "swaps": [
            # D30 (GAT) -> N30 (AAT) inside "Asp-Thr-Gly" active site
            ("GATACAGG", "AATACAGG"),
            ("GATACGGG", "AATACGGG"),
        ]
    }
]

def load_dataset(filepath):
    data = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            data.append(json.loads(line))
    return data

def extract_dna(json_obj):
    """Extracts raw DNA from the user prompt."""
    content = json_obj['messages'][0]['content']
    match = re.search(r"<dna_start>(.*?)<dna_end>", content)
    return match.group(1) if match else None

def inject_mutation(dna_seq, mutation_def):
    """
    Tries multiple wildcard patterns to find a match.
    Returns: (Success_Bool, New_DNA_Seq)
    """
    for wt, mut in mutation_def["swaps"]:
        if wt in dna_seq:
            # Perform the swap (only the first occurrence to avoid destroying structure)
            new_seq = dna_seq.replace(wt, mut, 1)
            return True, new_seq
            
    return False, dna_seq

def update_reasoning_and_report(original_msg, mutation_def):
    """Rewrites the Assistant's response to reflect the new mutation."""
    assistant_content = original_msg['messages'][1]['content']
    
    # 1. Parse existing Think vs Report
    think_match = re.search(r"<think>(.*?)</think>", assistant_content, re.DOTALL)
    if not think_match: return None
    
    original_think = think_match.group(1)
    original_report = assistant_content.split("</think>")[-1]

    # 2. Modify Reasoning Trace
    new_step_2 = f"Step 2: Identified surveillance mutations: {mutation_def['name']}."
    new_step_4 = f"Step 4: Significant resistance penalties detected for: {', '.join(mutation_def['drugs'])}."
    
    new_think = re.sub(r"Step 2:.*?\.", new_step_2, original_think)
    new_think = re.sub(r"Step 4:.*?\.", new_step_4, new_think)
    
    # 3. Modify Clinical Report
    new_report = original_report.replace(
        "No major surveillance mutations detected.", 
        f"{mutation_def['name']}"
    )
    
    # Update Drug Scores
    for drug in mutation_def['drugs']:
        # Escape parenthesis for regex: "LMV (3TC)" -> "LMV \(3TC\)"
        drug_esc = re.escape(drug)
        # Look for "Susceptible" entry and replace
        new_report = re.sub(
            f"- {drug_esc}: Susceptible", 
            f"- {drug}: {mutation_def['description']}", 
            new_report
        )

    # Reassemble
    new_content = f"<think>{new_think}</think>{new_report}"
    
    new_msg = original_msg.copy()
    new_msg['messages'] = [
        original_msg['messages'][0], 
        {"role": "assistant", "content": new_content}
    ]
    return new_msg

def main():
    print(f"Loading {INPUT_FILE}...")
    dataset = load_dataset(INPUT_FILE)
    
    # 1. Scan Current Stats
    susceptible_indices = []
    resistant_count = 0
    
    for i, entry in enumerate(dataset):
        report = entry['messages'][1]['content']
        # Simple heuristic to classify existing data
        if "High-Level" in report or "Intermediate" in report:
            resistant_count += 1
        else:
            susceptible_indices.append(i)
            
    print(f"Original Stats: {resistant_count} Resistant, {len(susceptible_indices)} Susceptible.")
    
    needed = TARGET_TOTAL_RESISTANT - resistant_count
    if needed <= 0:
        print("Dataset is already sufficiently resistant.")
        return

    print(f"Goal: Generate {needed} new resistant examples.")
    
    # 2. Augmentation Loop
    augmented_data = []
    random.shuffle(susceptible_indices)
    
    processed_count = 0
    idx_ptr = 0
    
    # Progress Bar
    pbar = tqdm(total=needed, desc="Injecting Mutations")
    
    while processed_count < needed and idx_ptr < len(susceptible_indices):
        idx = susceptible_indices[idx_ptr]
        original_entry = dataset[idx]
        dna = extract_dna(original_entry)
        
        if not dna: 
            idx_ptr += 1
            continue
            
        # Try mutations in random order until one works
        random.shuffle(KNOWN_MUTATIONS)
        mutation_success = False
        
        for mutation in KNOWN_MUTATIONS:
            success, new_dna = inject_mutation(dna, mutation)
            if success:
                new_entry = update_reasoning_and_report(original_entry, mutation)
                if new_entry:
                    # Update DNA in User Prompt
                    old_prompt = new_entry['messages'][0]['content']
                    new_prompt = re.sub(r"<dna_start>.*?<dna_end>", f"<dna_start>{new_dna}<dna_end>", old_prompt)
                    new_entry['messages'][0]['content'] = new_prompt
                    
                    augmented_data.append(new_entry)
                    processed_count += 1
                    pbar.update(1)
                    mutation_success = True
                    break # Stop trying mutations for this sequence
        
        # If we failed to inject ANY mutation, we just move to the next sequence
        # (This is expected for some weird subtypes)
        idx_ptr += 1
        
    pbar.close()
    
    # 3. Save
    final_dataset = dataset + augmented_data
    random.shuffle(final_dataset)
    
    print(f"\nSuccess: Augmented {len(augmented_data)} examples.")
    print(f"Saving total {len(final_dataset)} examples to {OUTPUT_FILE}...")
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        for entry in final_dataset:
            f.write(json.dumps(entry) + "\n")

if __name__ == "__main__":
    main()