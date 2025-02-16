#!/usr/bin/env python3
"""
PanSSRAtor – A Pan‑Species SSR Annotator

Usage:
  Genome Mode:
    python main.py --mode genome --genome_dir ./genomes/ --annot_dir ./annotations/ --output markers.tsv
  Genotype Mode:
    python main.py --mode genotype --reference ref.fasta --markers markers.tsv --bam_dir ./bams/ --output genotypes.csv
"""

import argparse
import time
from panssrator import config, utils, io_tools, ssr_discovery, annotator, primer_design, epcr, genotyper, marker_filter

def genome_mode(genome_dir: str, annot_dir: str, output: str):
    utils.logger.info("Running Genome Mode")
    pairs = io_tools.get_genome_annotation_pairs(genome_dir, annot_dir)
    all_markers = []
    for genome_file, annot_file in pairs:
        utils.logger.info("Processing genome: %s", genome_file)
        # For each FASTA in the genome file (could be multiple contigs)
        for header, seq in io_tools.read_fasta(genome_file):
            # Detect SSRs – if available, you may add chromosome information here.
            ssrs = ssr_discovery.detect_ssrs(seq)
            # Load annotation and add chromosome field
            annot_trees = annotator.load_annotation(annot_file)
            for rec in ssrs:
                rec["chrom"] = header  # assume header equals chromosome ID
                rec["annotation"] = annotator.annotate_ssr(rec, annot_trees)
                # Design primers: extract flanking region from the sequence
                rec["primers"] = primer_design.design_primers_for_ssr(rec, seq, flank=config.FLANK_SIZE)
                # Run ePCR simulation
                if rec.get("primers"):
                    # Extract primer pair (using the first left and right primer, if available)
                    primer_pair = {
                        "forward": rec["primers"].get("PRIMER_LEFT_0_SEQUENCE", ""),
                        "reverse": rec["primers"].get("PRIMER_RIGHT_0_SEQUENCE", "")
                    }
                    rec["amplicon_sizes"] = epcr.simulate_epcr(seq, primer_pair, max_cost=config.MAX_EPCR_COST)
            all_markers.extend(ssrs)
    # Filter markers
    filtered_markers = marker_filter.filter_markers(all_markers)
    utils.logger.info("Total markers detected: %d; Filtered markers: %d", len(all_markers), len(filtered_markers))
    # Write results to output (as TSV)
    with open(output, "w") as f:
        headers = ["chrom", "start", "end", "motif", "repeat_count", "annotation", "primers", "amplicon_sizes"]
        f.write("\t".join(headers) + "\n")
        for rec in filtered_markers:
            line = [str(rec.get(col, "")) for col in headers]
            f.write("\t".join(line) + "\n")
    utils.logger.info("Marker database saved to %s", output)

def genotype_mode(reference: str, markers_file: str, bam_dir: str, output: str):
    utils.logger.info("Running Genotype Mode")
    # Load reference genome – assume a single FASTA file; you can extend to use an indexed version.
    ref_seqs = {}
    for header, seq in io_tools.read_fasta(reference):
        ref_seqs[header] = seq
    # Load marker file (TSV)
    markers = []
    with open(markers_file, "r") as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            marker = {
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "motif": parts[3],
                "repeat_count": int(parts[4])
                # Additional fields can be added as needed.
            }
            markers.append(marker)
    # List BAM files
    bam_files = io_tools.list_files_in_dir(bam_dir, extensions=[".bam"])
    genotype_calls = {}
    for bam in bam_files:
        utils.logger.info("Processing BAM file: %s", bam)
        genotype_calls[bam] = {}
        for marker in markers:
            # Use the reference sequence for the marker’s chromosome.
            seq = ref_seqs.get(marker["chrom"], "")
            if not seq:
                continue
            gt = genotyper.genotype_marker(bam, marker)
            genotype_calls[bam][marker["start"]] = gt
    # Write genotype table as CSV
    with open(output, "w") as f:
        header_line = "BAM_file,Marker_start,Genotype\n"
        f.write(header_line)
        for bam, calls in genotype_calls.items():
            for marker_start, call in calls.items():
                f.write(f"{bam},{marker_start},{call}\n")
    utils.logger.info("Genotype calls saved to %s", output)

def parse_args():
    parser = argparse.ArgumentParser(description="PanSSRAtor – Pan‑Species SSR Annotator")
    parser.add_argument("--mode", choices=["genome", "genotype"], required=True,
                        help="Select the mode of operation: genome (SSR discovery) or genotype (genotyping from BAM)")
    parser.add_argument("--genome_dir", help="Directory of genome FASTA files (for genome mode)")
    parser.add_argument("--annot_dir", help="Directory of annotation (GFF/GTF) files (for genome mode)")
    parser.add_argument("--reference", help="Reference genome FASTA (for genotype mode)")
    parser.add_argument("--markers", help="Marker file (from genome mode, for genotype mode)")
    parser.add_argument("--bam_dir", help="Directory of BAM files (for genotype mode)")
    parser.add_argument("--output", required=True, help="Output file (or prefix) for results")
    return parser.parse_args()

def main():
    args = parse_args()
    start = time.time()
    if args.mode == "genome":
        if not args.genome_dir or not args.annot_dir:
            utils.do_error("Genome mode requires --genome_dir and --annot_dir.")
        genome_mode(args.genome_dir, args.annot_dir, args.output)
    elif args.mode == "genotype":
        if not args.reference or not args.markers or not args.bam_dir:
            utils.do_error("Genotype mode requires --reference, --markers, and --bam_dir.")
        genotype_mode(args.reference, args.markers, args.bam_dir, args.output)
    end = time.time()
    utils.logger.info("PanSSRAtor run time: %.2f minutes", (end - start) / 60)

if __name__ == '__main__':
    main()

