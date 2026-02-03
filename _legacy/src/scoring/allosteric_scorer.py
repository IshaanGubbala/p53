from __future__ import annotations

import re

# Known allosteric hotspots in p53 Core Domain
ALLOSTERIC_SITES = {
    "Y220_POCKET": [220, 221, 222, 223, 229, 230],
    "L1_S3_POCKET": [112, 113, 114, 115, 116], 
    "DNA_BINDING": [248, 273, 280], 
}

def get_distance_to_site(pos: int, site_residues: list[int], dist_map: dict[tuple[int, int], float]) -> float:
    """Compute the minimum distance from a position to any residue in the site."""
    if pos in site_residues:
        return 0.0
    
    min_dist = float("inf")
    for res in site_residues:
        d = dist_map.get((pos, res))
        if d is None:
             d = dist_map.get((res, pos))
        if d is not None:
            min_dist = min(min_dist, d)
            
    return min_dist

def allosteric_distance_score(
    mutations: list[str], 
    dist_map: dict[tuple[int, int], float],
    sites: list[str] = ["Y220_POCKET"]
) -> float:
    """
    Score based on proximity to allosteric sites.
    Returns: Average distance to the nearest allosteric site.
    """
    if not mutations:
        return 0.0
        
    total_dist = 0.0
    count = 0
    
    target_residues = []
    for site_name in sites:
        target_residues.extend(ALLOSTERIC_SITES.get(site_name, []))
        
    for mut in mutations:
        match = re.search(r"(\d+)", mut)
        if not match:
            continue
        pos = int(match.group(1))
        
        d = get_distance_to_site(pos, target_residues, dist_map)
        if d != float("inf"):
            total_dist += d
            count += 1
            
    if count == 0:
        return 999.0 
        
    return total_dist / count
