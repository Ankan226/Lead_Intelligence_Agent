from src.constants import KEYWORDS, HUB_LOCATIONS, ROLE_KEYWORDS, CURRENT_YEAR

def score_lead(row):
    score = 0

    title = str(row["title"]).lower()
    if sum(k in title for k in KEYWORDS) >= 2:
        score += 40

    if row["year"]:
        age = CURRENT_YEAR - row["year"]
        if age <= 2:
            score += 40
        elif age <= 5:
            score += 20

    affiliation = str(row.get("affiliation", "")).lower()
    if any(r in affiliation for r in ROLE_KEYWORDS):
        score += 30

    location = str(row.get("person_location", "")).lower()
    if any(h in location for h in HUB_LOCATIONS):
        score += 10

    return min(score, 100)
