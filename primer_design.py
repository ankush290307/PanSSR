# panssrator/primer_design.py
import primer3
from panssrator import config, utils

def design_primers_for_ssr(ssr_record: dict, genome_seq: str, flank: int = config.FLANK_SIZE,
                           custom_params: dict = None) -> dict:
    """
    Extract flanking sequence around an SSR locus and design primers using Primer3.
    
    Parameters:
      ssr_record: Dictionary with SSR information (must include 'start' and 'end').
      genome_seq: The entire sequence (e.g., contig or chromosome) in which the SSR is located.
      flank: Number of bases to extract upstream and downstream.
      custom_params: Optional dictionary of primer design parameters (overrides defaults).
    
    Returns:
      A dictionary of primer design results.
    """
    start = max(0, ssr_record["start"] - flank - 1)  # convert to 0-indexed
    end = min(len(genome_seq), ssr_record["end"] + flank)
    template_seq = genome_seq[start:end]
    
    # Define the target region within the template
    target_start = flank  # because we extracted 'flank' bases upstream
    target_length = ssr_record["end"] - ssr_record["start"] + 1
    
    seq_args = {
        'SEQUENCE_ID': f"SSR_{ssr_record['start']}_{ssr_record['end']}",
        'SEQUENCE_TEMPLATE': template_seq,
        'SEQUENCE_TARGET': [target_start, target_length]
    }
    
    params = custom_params if custom_params else config.PRIMER_PARAMS
    
    try:
        result = primer3.designPrimers(seq_args, params)
    except Exception as e:
        utils.logger.error("Primer3 design failed: %s", e)
        result = {}
    return result

if __name__ == '__main__':
    # Example: design primers for a dummy SSR in a synthetic genome sequence.
    dummy_seq = "N" * config.FLANK_SIZE + "AT" * 10 + "N" * config.FLANK_SIZE
    ssr = {"start": config.FLANK_SIZE + 1, "end": config.FLANK_SIZE + 20, "motif": "AT"}
    primers = design_primers_for_ssr(ssr, dummy_seq)
    utils.logger.info("Primer design results: %s", primers)

