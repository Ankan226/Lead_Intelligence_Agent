import os
import streamlit as st
import pandas as pd

from src.pubmed_client import fetch_pubmed_leads
from src.enrichment import enrich_lead
from src.scoring import score_lead

os.makedirs("output", exist_ok=True)

st.set_page_config(
    page_title="3D In-Vitro Lead Qualification Dashboard",
    layout="wide"
)

st.title("3D In-Vitro Lead Qualification Dashboard")
st.caption(
    "Identifies, enriches, and ranks life-science professionals "
    "based on their probability of working with 3D in-vitro models."
)



st.sidebar.header(" Filters")
search = st.sidebar.text_input("Search (Name, Title, Company, Location)")
min_score = st.sidebar.slider("Minimum Probability Score", 0, 100, 0, 5)

query = st.text_input("PubMed search keywords", "liver toxicity 3D in vitro")

if st.button("Generate Leads"):
    with st.spinner("Processing leads..."):
        df = fetch_pubmed_leads(query)

        enriched = []
        for _, row in df.iterrows():
            enriched.append({**row.to_dict(), **enrich_lead(row)})

        df = pd.DataFrame(enriched)
        df["probability_score"] = df.apply(score_lead, axis=1)
        df = df.sort_values("probability_score", ascending=False)
        df["rank"] = range(1, len(df) + 1)

        if search:
            q = search.lower()
            df = df[
                df["name"].str.lower().str.contains(q)
                | df["title"].str.lower().str.contains(q)
                | df["company"].str.lower().str.contains(q)
                | df["person_location"].str.lower().str.contains(q)
            ]

        df = df[df["probability_score"] >= min_score]

        final = df[
            [
                "rank",
                "probability_score",
                "name",
                "title",
                "company",
                "person_location",
                "company_hq",
                "work_mode",
                "email"
            ]
        ]

        final.columns = [
            "Rank",
            "Probability Score",
            "Name",
            "Title",
            "Company",
            "Person Location",
            "Company HQ",
            "Work Mode",
            "Email"
        ]

        st.dataframe(final, use_container_width=True)

        final.to_csv("output/ranked_leads.csv", index=False)

        st.download_button(
            "Download Qualified Leads (CSV)",
            data=final.to_csv(index=False),
            file_name="ranked_leads.csv",
            mime="text/csv"
        )
