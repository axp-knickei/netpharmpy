import requests
import json
import time

def test_endpoint(name, url):
    print(f"\n{'='*60}")
    print(f"TESTING: {name}")
    print(f"URL:     {url}")
    print(f"{'-'*60}")
    
    try:
        start_time = time.time()
        response = requests.get(url, timeout=15, headers={'Accept': 'application/json'})
        elapsed = time.time() - start_time
        
        print(f"Status:  {response.status_code}")
        print(f"Time:    {elapsed:.2f}s")
        
        if response.status_code == 200:
            data = response.json()
            
            # Analyze structure
            if isinstance(data, list):
                print(f"Result:  Found {len(data)} items")
                if len(data) > 0:
                    print("\n[First Item Sample]:")
                    print(json.dumps(data[0], indent=2)[:500] + "...")
                    
                    # Check for UniProt/Gene info
                    has_uniprot = any(d.get('databaseName') == 'UniProt' for d in data)
                    print(f"\nContains UniProt Entities? {'YES' if has_uniprot else 'NO'}")
            else:
                print("Result:  Returned object (not a list)")
                print(json.dumps(data, indent=2)[:500])
        else:
            print("FAILED. Response text:")
            print(response.text[:500])
            
    except Exception as e:
        print(f"CRITICAL ERROR: {str(e)}")

# Test Case: "Cytokine Signaling in Immune system"
PATHWAY_ID = "R-HSA-1280215"

print(f"DIAGNOSING PATHWAY: {PATHWAY_ID}")

# 1. The Method Currently Failing (likely returns Complexes, not Proteins)
test_endpoint(
    "Current Method (Participants)", 
    f"https://reactome.org/ContentService/data/participants/{PATHWAY_ID}"
)

# 2. The Proposed Fix (Reference Entities - Should return flat list)
test_endpoint(
    "Proposed Fix (Participating Reference Entities)", 
    f"https://reactome.org/ContentService/data/event/{PATHWAY_ID}/participatingReferenceEntities"
)

# 3. Alternative Method (Physical Entities)
test_endpoint(
    "Alternative Method (Participating Physical Entities)", 
    f"https://reactome.org/ContentService/data/event/{PATHWAY_ID}/participatingPhysicalEntities"
)