import os
import streamlit as st
import pandas as pd

from src.pubmed_client import fetch_pubmed_leads
from src.enrichment import enrich_lead
from src.scoring import score_lead

# ✅ AUTO-CREATE OUTPUT FOLDER
os.makedirs("output", exist_ok=True)

st.set_page_config(page_title="Lead Intelligence Agent", layout="wide")
st.title(" Scientific Lead Intelligence Dashboard")

query = st.text_input(
    "Enter PubMed search keywords",
    "liver toxicity 3D in vitro"
)

if st.button("Generate & Rank Leads"):
    with st.spinner("Fetching live data from PubMed..."):
        df = fetch_pubmed_leads(query)

        if df.empty:
            st.warning("No leads found.")
        else:
            enriched_rows = []
            for _, row in df.iterrows():
                enriched_rows.append({**row.to_dict(), **enrich_lead(row)})

            df = pd.DataFrame(enriched_rows)

            df["Probability Score"] = df.apply(score_lead, axis=1)
            df = df.sort_values("Probability Score", ascending=False)
            df["Rank"] = range(1, len(df) + 1)

            df = df[
                [
                    "Rank",
                    "name",
                    "paper_title",
                    "year",
                    "company_type",
                    "person_location",
                    "company_hq",
                    "email",
                    "Probability Score"
                ]
            ]

            df.columns = [
                "Rank",
                "Name",
                "Paper Title",
                "Year",
                "Company Type",
                "Person Location",
                "Company HQ",
                "Email",
                "Probability Score"
            ]

            st.dataframe(df, use_container_width=True)

            df.to_csv("output/ranked_leads.csv", index=False)

            st.download_button(
                "Download CSV",
                data=df.to_csv(index=False),
                file_name="ranked_leads.csv",
                mime="text/csv"
            )
