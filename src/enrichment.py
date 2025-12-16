def enrich_lead(row):
    aff = str(row.get("affiliation", "")).lower()

    if "university" in aff or "institute" in aff:
        company = "Academic"
    elif any(x in aff for x in ["pharma", "biotech", "inc", "ltd"]):
        company = "Pharma/Biotech"
    else:
        company = "Research Org"

    location = "Unknown"
    for city in ["boston", "cambridge", "basel", "london", "san francisco"]:
        if city in aff:
            location = city.title()
            break

    email = None
    name_parts = row["name"].split()
    if len(name_parts) >= 2 and company != "Academic":
        email = f"{name_parts[0].lower()}.{name_parts[-1].lower()}@company.com"

    return {
        "company": company,
        "person_location": location,
        "company_hq": location,
        "work_mode": "Remote" if "remote" in aff else "Onsite",
        "email": email
    }
