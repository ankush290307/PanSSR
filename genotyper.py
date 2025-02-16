# panssrator/genotyper.py
import pysam
import re
from collections import Counter
from panssrator import config, utils

def count_repeat_units(seq: str, motif: str) -> int:
    """
    Count the number of times a given motif is repeated consecutively in a sequence.
    This simple approach uses regex matching.
    """
    # Build a regex pattern for the motif repeated consecutively
    pattern = f"(?:{motif})+"
    match = re.search(pattern, seq, re.IGNORECASE)
    if match:
        repeated_seq = match.group(0)
        return len(repeated_seq) // len(motif)
    return 0

def genotype_marker(bam_file: str, ssr_record: dict) -> dict:
    """
    Extract reads from a BAM file that overlap the SSR region,
    count the repeat units in each read, and call a genotype.
    
    Parameters:
      bam_file: Path to the BAM file.
      ssr_record: Dictionary with keys including 'chrom', 'start', 'end', 'motif'.
    
    Returns:
      A dictionary containing the most common repeat count and supporting read counts.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    allele_counts = Counter()
    chrom = ssr_record.get("chrom", None)
    if not chrom:
        utils.do_error("SSR record does not contain chromosome information.")
    start = ssr_record["start"]
    end = ssr_record["end"]
    for read in bam.fetch(chrom, start - 1, end + 1):
        if read.mapping_quality < config.MIN_MAPQ:
            continue
        # Determine the alignment region overlapping the SSR.
        # This is simplified: in practice, you may want to extract the portion
        # of the read corresponding to the SSR region from the read’s CIGAR string.
        read_seq = read.query_sequence
        repeat_count = count_repeat_units(read_seq, ssr_record["motif"])
        if repeat_count:
            allele_counts[repeat_count] += 1
    bam.close()
    
    # Determine genotype based on allele counts
    total = sum(allele_counts.values())
    if total < config.MIN_READ_SUPPORT:
        genotype = None  # Insufficient read support
    else:
        # For diploid species, expect at most two alleles.
        most_common = allele_counts.most_common(2)
        if len(most_common) == 1:
            genotype = {"alleles": [most_common[0][0]], "support": most_common[0][1]}
        else:
            genotype = {"alleles": [most_common[0][0], most_common[1][0]],
                        "support": {most_common[0][0]: most_common[0][1],
                                    most_common[1][0]: most_common[1][1]}}
    return {"allele_counts": dict(allele_counts), "genotype": genotype}

if __name__ == '__main__':
    # Example test: Replace 'example.bam' with an actual BAM file to run this test.
    ssr_example = {"chrom": "chr1", "start": 100, "end": 140, "motif": "AT"}
    gt = genotype_marker("example.bam", ssr_example)
    utils.logger.info("Genotype call: %s", gt)

