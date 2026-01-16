import requests
import json

# 1. The Endpoints to test (Main and Backup)
URLS = [
    "https://hivdb.stanford.edu/graphql",
    "https://hivdb.stanford.edu/sierra/graphql" 
]

# 2. A known 'Gold Standard' short sequence (HIV-1 Protease/RT start)
# This removes your data from the equation to test the connection first.
TEST_SEQ = "CCTCAAATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTT"

# 3. The Query (Standard Sierra)
QUERY = """
query ($sequences: [SequenceInput]!) {
  sequenceAnalysis(sequences: $sequences) {
    sequence { header }
    drugResistance {
      mutationsByGene {
        gene { name }
        mutations { text }
      }
    }
  }
}
"""

def test_connection():
    payload = {
        "query": QUERY,
        "variables": {
            "sequences": [
                {"header": "TEST_SEQ_001", "sequence": TEST_SEQ}
            ]
        }
    }

    print(f"--- Starting Diagnostic ---")
    
    for url in URLS:
        print(f"\nTesting URL: {url}")
        try:
            r = requests.post(url, json=payload, timeout=10)
            
            print(f"Status Code: {r.status_code}")
            
            if r.status_code == 200:
                print("SUCCESS! API Response:")
                print(json.dumps(r.json(), indent=2))
                return url # Return the working URL
            
            elif r.status_code == 400:
                print("FAIL (400). Server says:")
                # THIS IS THE KEY: The server usually tells you WHY it failed here
                print(r.text)
                
            elif r.status_code == 405:
                print("FAIL (405). Method Not Allowed (Did you use GET instead of POST?)")
                
            else:
                print(f"FAIL ({r.status_code}). Response: {r.text}")

        except Exception as e:
            print(f"Connection Error: {e}")

    return None

if __name__ == "__main__":
    working_url = test_connection()
    if working_url:
        print(f"\n\n>>> RECOMMENDATION: Update CONFIG['API_URL'] to: {working_url}")
    else:
        print("\n\n>>> CONCLUSION: Both endpoints failed. Check the error messages above.")