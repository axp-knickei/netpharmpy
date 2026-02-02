"""
Manual test script for DAVID Web Service integration.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from netpharm.enrichment.david import DavidEnrichmentClient

def test_david_connection(email):
    print(f"--- Testing DAVID Connection for: {email} ---")
    client = DavidEnrichmentClient(email)
    
    try:
        print("1. Authenticating...")
        is_auth = client.authenticate()
        print(f"   Success: {is_auth}")
        
        print("\n2. Submitting small gene list (TP53, TNF)...")
        # Human Entrez IDs for TP53 (7157) and TNF (7124)
        genes = ["7157", "7124"]
        response = client.add_gene_list(genes, list_name="Test_NetPharmPy", list_type="ENTREZ_GENE_ID")
        print(f"   List Submission Response: {response[:100]}...") # Print first 100 chars
        
        print("\n3. Setting species to Homo sapiens (9606)...")
        client.set_categories(["GOTERM_BP_DIRECT"])
        
        print("\n4. Requesting Chart Report...")
        results = client.get_chart_report(threshold=1.0, min_count=1) # High threshold to ensure we get something
        
        print(f"\n--- RESULTS FOUND: {len(results)} ---")
        for res in results[:5]:
            print(f"   [{res.get('Category')}] {res.get('Term')}: P={res.get('PValue')}")
            
        print("\n✓ DAVID Integration Test Passed!")
        
    except Exception as e:
        print(f"\n❌ DAVID Integration Test Failed!")
        print(f"   Error: {str(e)}")
        print("\n💡 Possible reasons:")
        print("   - Email is not registered with DAVID Web Service.")
        print("   - Network/Firewall is blocking SOAP requests to david.ncifcrf.gov.")
        print("   - DAVID server is currently down or under maintenance.")

if __name__ == "__main__":
    email = "alexprima@gmail.com"
    test_david_connection(email)
