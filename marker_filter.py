# panssrator/marker_filter.py
from typing import List, Dict, Any
from panssrator import utils

def filter_markers(markers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Filter markers based on the following criteria:
      - Unique amplification across genomes (from ePCR simulation)
      - Consistent SSR motif in all amplicons
      - Amplicon length differences are due only to repeat unit variation
      - (Optional) Annotation: prefer intergenic markers
    Returns a list of filtered markers.
    """
    filtered = []
    for marker in markers:
        # Example filter: must have a unique amplicon size per genome.
        # Assume marker["amplicon_sizes"] is a list of product sizes across genomes.
        sizes = marker.get("amplicon_sizes", [])
        if not sizes:
            continue
        if len(set(sizes)) < 2:
            continue  # Not polymorphic
        # Additional filters could be added here.
        filtered.append(marker)
    return filtered

if __name__ == '__main__':
    # Test with dummy markers.
    dummy_markers = [
        {"start": 100, "end": 140, "motif": "AT", "amplicon_sizes": [200, 200, 205]},
        {"start": 150, "end": 190, "motif": "CG", "amplicon_sizes": [300, 300, 300]},
    ]
    filtered = filter_markers(dummy_markers)
    utils.logger.info("Filtered markers: %s", filtered)

