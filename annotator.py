# panssrator/annotator.py
import os
from typing import List, Dict, Optional
from intervaltree import Interval, IntervalTree
from panssrator import utils

def load_annotation(annotation_file: str) -> Dict[str, IntervalTree]:
    """
    Load a GFF/GTF file and return a dictionary mapping chromosomes to an IntervalTree of features.
    
    Each feature is stored as a dictionary with keys: 'type', 'start', 'end', 'strand', 'attributes'.
    """
    trees = {}
    with open(annotation_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            feature = {
                "type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "attributes": attributes
            }
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom][start:end+1] = feature  # end+1 because intervaltree is half-open [start, end)
    return trees

def annotate_ssr(ssr_record: dict, annot_trees: Dict[str, IntervalTree]) -> Optional[dict]:
    """
    Given an SSR record and annotation interval trees, find the first overlapping feature.
    Returns the feature dictionary if found; otherwise, None.
    """
    chrom = ssr_record.get("chrom", None)
    if not chrom or chrom not in annot_trees:
        return None
    overlaps = annot_trees[chrom].overlap(ssr_record["start"], ssr_record["end"]+1)
    if overlaps:
        # Return the first overlapping feature (or apply a priority rule if needed)
        return list(overlaps)[0].data
    return None

if __name__ == '__main__':
    # Example: assume an annotation file "example.gff" exists.
    try:
        trees = load_annotation("example.gff")
        test_ssr = {"chrom": "chr1", "start": 100, "end": 140, "motif": "AT"}
        feature = annotate_ssr(test_ssr, trees)
        if feature:
            utils.logger.info("SSR is located in feature: %s", feature)
        else:
            utils.logger.info("SSR not located in any annotated feature.")
    except Exception as e:
        utils.do_error(str(e))

