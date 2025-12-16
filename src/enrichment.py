def enrich_lead(row):
    affiliation = str(row.get("affiliation", "")).lower()

    # Company type inference
    if "university" in affiliation or "institute" in affiliation:
        company_type = "Academic"
    elif any(x in affiliation for x in ["pharma", "inc", "ltd", "corp"]):
        company_type = "Pharma/Biotech"
    else:
        company_type = "Research Organization"

    # Location inference
    person_location = "Unknown"
    for city in ["boston", "cambridge", "basel", "london", "san francisco"]:
        if city in affiliation:
            person_location = city.title()
            break

    company_hq = person_location

    # Email heuristic (explicitly a demo assumption)
    email = None
    if company_type != "Academic":
        parts = row["name"].split()
        if len(parts) >= 2:
            email = f"{parts[0].lower()}.{parts[-1].lower()}@company.com"

    return {
        "company_type": company_type,
        "person_location": person_location,
        "company_hq": company_hq,
        "email": email
    }
