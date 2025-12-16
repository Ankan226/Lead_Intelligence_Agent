from src.constants import KEYWORDS, HUB_LOCATIONS, ROLE_KEYWORDS, CURRENT_YEAR

def score_lead(row):
    score = 0

    # Scientific intent
    title = row["paper_title"].lower()
    keyword_hits = sum(1 for k in KEYWORDS if k in title)
    if keyword_hits >= 2:
        score += 40
    elif keyword_hits == 1:
        score += 20

    # Recency
    if row["year"]:
        age = CURRENT_YEAR - row["year"]
        if age <= 2:
            score += 40
        elif age <= 5:
            score += 20

    # Role fit
    affiliation = str(row.get("affiliation", "")).lower()
    if any(role in affiliation for role in ROLE_KEYWORDS):
        score += 30

    # Location hub
    location = str(row.get("person_location", "")).lower()
    if any(hub in location for hub in HUB_LOCATIONS):
        score += 10

    return min(score, 100)
