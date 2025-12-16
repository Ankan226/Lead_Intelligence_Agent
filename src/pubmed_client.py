from Bio import Entrez
import pandas as pd

Entrez.email = "ankanpal255@gmail.com"

def fetch_pubmed_leads(query, max_results=20):
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_results,
        sort="relevance"
    )
    record = Entrez.read(handle)
    ids = record["IdList"]

    leads = []

    for pid in ids:
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=pid,
            rettype="medline",
            retmode="text"
        )
        data = fetch_handle.read()

        title = ""
        year = None
        affiliation = ""
        authors = []
        corresponding = False

        for line in data.split("\n"):
            if line.startswith("TI  -"):
                title = line.replace("TI  -", "").strip()

            if line.startswith("DP  -"):
                try:
                    year = int(line[6:10])
                except:
                    year = None

            if line.startswith("FAU -"):
                authors.append(line.replace("FAU -", "").strip())

            if line.startswith("AD  -"):
                affiliation = line.replace("AD  -", "").strip()

            if "correspond" in line.lower():
                corresponding = True

        for author in authors:
            leads.append({
                "name": author,
                "title": title,
                "year": year,
                "affiliation": affiliation,
                "corresponding_author": corresponding
            })

    return pd.DataFrame(leads)
