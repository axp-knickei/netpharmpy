"""
DAVID Web Service Client (SOAP via Requests)

Implements interaction with DAVID bioinformatics resources 
using standard SOAP/XML protocols over HTTP.
"""

import requests
import xml.etree.ElementTree as ET
import logging
import time
import re

class DavidEnrichmentClient:
    """
    Client for DAVID Web Service (SOAP).
    """
    
    WSDL_URL = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService"
    NAMESPACE = "http://service.session.sample"
    
    def __init__(self, email):
        self.email = email
        self.headers = {
            'Content-Type': 'text/xml; charset=utf-8'
        }
        self.session = requests.Session()
        
    def _send_soap(self, method, params=None):
        """
        Send raw SOAP request.
        """
        body_content = ""
        if params:
            for key, value in params.items():
                body_content += f"<{key}>{value}</{key}>"
                
        envelope = f"""<soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:ser="{self.NAMESPACE}">
   <soapenv:Header/>
   <soapenv:Body>
      <ser:{method}>
         {body_content}
      </ser:{method}>
   </soapenv:Body>
</soapenv:Envelope>"""

        response = self.session.post(self.WSDL_URL, data=envelope, headers=self.headers)
        
        if response.status_code != 200:
            raise Exception(f"SOAP Request Failed: {response.status_code} - {response.text}")
            
        return response.text

    def authenticate(self):
        """Authenticate with email."""
        response = self._send_soap("authenticate", {"arg0": self.email})
        if "true" not in response.lower():
            raise Exception("Authentication failed. Check email registration with DAVID.")
        return True

    def add_gene_list(self, genes, list_name="NetPharmPy_List", list_type="ENTREZ_GENE_ID"):
        """Add gene list to session."""
        gene_str = ",".join(map(str, genes))
        params = {
            "arg0": gene_str,
            "arg1": list_type,
            "arg2": list_name,
            "arg3": 0  # 0 = Gene list
        }
        response = self._send_soap("addList", params)
        # Parse return value (overlap count)
        # Assuming simple check for now
        return response

    def set_categories(self, categories):
        """Set categories (comma separated string)."""
        cat_str = ",".join(categories)
        self._send_soap("setCategories", {"arg0": cat_str})

    def get_chart_report(self, threshold=0.05, min_count=2):
        """Get chart report records."""
        params = {
            "arg0": threshold,
            "arg1": min_count
        }
        xml_response = self._send_soap("getChartReport", params)
        return self._parse_chart_report(xml_response)

    def _parse_chart_report(self, xml_text):
        """Parse the complex object return from getChartReport."""
        # Simple parser for DAVID SOAP response
        records = []
        try:
            # Remove namespaces to simplify finding tags
            xml_text = re.sub(r'xmlns:.*?"[^"]+"', '', xml_text)
            xml_text = re.sub(r'\w+:', '', xml_text) # remove prefix:
            
            root = ET.fromstring(xml_text)
            
            # Look for all items that have termName
            # In SOAP responses, these are usually <item> tags inside the return array
            
            # Strategy: Find all elements that have a 'termName' child
            # then extract siblings
            
            for item in root.findall(".//termName/.."):
                record = {}
                record['Category'] = item.findtext("categoryName")
                record['Term'] = item.findtext("termName")
                record['Count'] = item.findtext("listHits")
                record['%'] = item.findtext("percent")
                record['PValue'] = item.findtext("ease")
                record['Benjamini'] = item.findtext("benjamini")
                record['FDR'] = item.findtext("fdr")
                
                records.append(record)
                
        except Exception as e:
            logging.error(f"Error parsing DAVID XML: {e}")
            
        return records
