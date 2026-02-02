# DAVID Enrichment Workflow

This document outlines the workflow and assumptions for using the DAVID Web Service integration in NetpharmPy.

## Workflow
1.  **Authentication**: The user provides a registered DAVID email address.
2.  **Gene List Upload**: The pipeline submits the filtered target list (e.g., from network hubs) to DAVID via SOAP.
3.  **Analysis**: The pipeline requests a Chart Report for specified categories (GO:BP, KEGG, Reactome).
4.  **Download**: Results are parsed and saved as CSV/TSV files.

## Constraints & Limitations
-   **Registration Required**: Users must register their email with the DAVID Web Service (https://david.ncifcrf.gov/webservice/register.htm) before running this pipeline.
-   **Rate Limits**: DAVID imposes strict limits (approx. 200 jobs/day/user). Do not run this inside a loop.
-   **No Discovery**: This module is intended for **confirmatory** analysis of small, high-confidence gene lists, not for exploratory scanning of thousands of sets.
-   **SOAP Protocol**: This implementation uses a custom SOAP client (`netpharm.enrichment.david`) to avoid `zeep`/`suds` dependencies. It assumes the DAVID WSDL structure remains stable.

## Reproducibility
-   DAVID database versions update periodically. The metadata file saved with results records the date of analysis.
-   Exact p-values may vary slightly between DAVID versions.
-   Ensure you use the same email/account for consistent access logs.