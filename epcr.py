# panssrator/epcr.py
import re
import sys
from panssrator import config, utils

try:
    import tre
except ImportError:
    utils.do_error("The 'tre' module is required for ePCR simulation. Please install it from https://github.com/laurikari/tre/")

def compile_primer_pattern(primer_seq: str) -> object:
    """
    Compile a primer sequence (with ambiguity codes replaced) into a TRE regex pattern.
    """
    regex_seq = replace_ambiguity_codes(primer_seq)
    try:
        pattern = tre.compile(regex_seq, tre.EXTENDED)
    except Exception as e:
        utils.do_error(f"Error compiling primer pattern: {e}")
    return pattern

def simulate_epcr(genome_seq: str, primer_pair: dict, max_cost: int = config.MAX_EPCR_COST) -> list:
    """
    Simulate in silico PCR by searching for primer binding sites in genome_seq.
    
    Parameters:
      genome_seq: The target genome sequence.
      primer_pair: Dictionary with keys 'forward' and 'reverse' containing primer sequences.
      max_cost: Maximum allowed fuzzy matching cost.
      
    Returns:
      A list of predicted amplicon sizes (in bp) where both primers are found in correct orientation.
    """
    # Compile both primer patterns using the TRE module
    forward_pat = tre.compile(replace_ambiguity_codes(primer_pair.get("forward", "")), tre.EXTENDED)
    reverse_pat = tre.compile(replace_ambiguity_codes(primer_pair.get("reverse", "")), tre.EXTENDED)
    # For the reverse primer, in a full implementation we’d search for its reverse complement.
    # Here we assume the provided reverse primer is already in the proper orientation.
    
    f_matches = forward_pat.findall(genome_seq, fuzzyness=max_cost)
    r_matches = reverse_pat.findall(genome_seq, fuzzyness=max_cost)
    
    # To get positions, we can use search() instead of findall()
    f_positions = []
    pos = 0
    while True:
        m = forward_pat.search(genome_seq, pos, fuzzyness=max_cost)
        if not m:
            break
        f_positions.append(m.start())
        pos = m.start() + 1
    r_positions = []
    pos = 0
    while True:
        m = reverse_pat.search(genome_seq, pos, fuzzyness=max_cost)
        if not m:
            break
        r_positions.append(m.start())
        pos = m.start() + 1

    amplicon_sizes = []
    # For each forward match, find a reverse match that is downstream.
    for f in f_positions:
        for r in r_positions:
            if r > f:
                size = r - f + len(primer_pair.get("reverse", ""))
                if size <= config.MAX_EPCR_PRODUCT:
                    amplicon_sizes.append(size)
    return amplicon_sizes

def replace_ambiguity_codes(seq: str) -> str:
    """Wrapper to call utils.replace_ambiguity_codes."""
    return utils.replace_ambiguity_codes(seq)

if __name__ == '__main__':
    # Example test for ePCR simulation
    test_genome = "N" * 100 + "ATGCGT" + "N" * 50 + "CATGCA" + "N" * 100
    primers = {"forward": "ATGCGT", "reverse": "CATGCA"}
    products = simulate_epcr(test_genome, primers)
    utils.logger.info("Predicted amplicon sizes: %s", products)

