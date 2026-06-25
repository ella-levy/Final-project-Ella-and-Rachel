import os
import math
from collections import Counter


FASTA_PATH = "data/antibodies_expanded.fasta"   # where to read sequences from
RESULTS_DIR = "results"  # where to save output files
STD_MULTIPLIER = 1.0  # threshold = mean + this * std
                          # 1.0 catches ~16% of positions (flexible)
                          # 2.0 catches ~3%  of positions (strict)
MIN_REGION_LEN  = 3 # a CDR region must be at least this many positions long

def parse_fasta(filepath):
    """
    Read FASTA file that contains multiple homo sepiens sequences.
    Returns a dict: { header_string : sequence_string }
    """

    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous sequence before starting a new one
                if current_header and current_seq:
                    sequences[current_header] = "".join(current_seq).upper()
                current_header = line[1:]   # remove the ">"
                current_seq = []
            else:
                # Sequence lines can be multiple lines
                current_seq.append(line)

        # Save the last sequence in the file
        if current_header and current_seq:
            sequences[current_header] = "".join(current_seq).upper()

    # print summary of what we loaded
    print(f"Loaded sequences from {filepath}:")
    for header, seq in sequences.items():
        parts = header.split("|")
        pdb_id = parts[0] if parts else header
        chain_desc = parts[2] if len(parts) > 2 else "?"
        print(f"{pdb_id:10s} | {chain_desc[:40]:40s} | {len(seq)} aa")
    return sequences

def calculate_entropy(sequences):
    """
    For each position (column) calculate the entropy
    The formula:  H = -sum( pᵢ * log2(pᵢ) )
    Input:  list of sequence strings
    Output: list of entropy values, one per position
    """

    if not sequences:
        return []

    # Use the shortest sequence so all sequences are comparable column by column
    num_positions = min(len(s) for s in sequences)
    num_sequences = len(sequences)
    entropy_values = []

    for col in range(num_positions):
        # Grab the amino acid from every sequence at this column
        column = [seq[col] for seq in sequences]

        # Count how many times each amino acid appears
        counts = Counter(column)

        # Calculate entropy for this column
        H = 0.0
        for aa, count in counts.items():
            if aa in ("-", "X"):
                continue   # skip gap for unknown characters

            p = count / num_sequences  # frequency of this amino acid here
            H -= p * math.log2(p)

        entropy_values.append(round(H, 5))

    return entropy_values


def find_cdr_regions(entropy_values, std_multiplier, min_region_length):
    """
    Find positions above the threshold then group consecutive ones into CDR regions
    Threshold = mean + std_multiplier * std

    Input:
        entropy_values - list of entropy per position
        std_multiplier - how strict the threshold is (1.0 = flexible, 2.0 = strict)
        min_region_length - minimum consecutive positions to count as a region

    Output: dict with threshold, positions above it, and grouped regions
    """

    n = len(entropy_values)
    mean = sum(entropy_values) / n
    variance = sum((x - mean) ** 2 for x in entropy_values) / n
    std = math.sqrt(variance)
    threshold = mean + std_multiplier * std

#print info
    print(f"\n Entropy statistics:")
    print(f"Mean: {mean:.4f}")
    print(f"Std: {std:.4f}")
    print(f"Threshold: {threshold:.4f}  (mean + {std_multiplier}*std)")

    # Find every position whose entropy is above the threshold
    high_positions = [i for i, h in enumerate(entropy_values) if h >= threshold]
    print(f"positions above threshold: {len(high_positions)} / {n}")

    # Group consecutive positions into regions
    regions = []
    if high_positions:
        start = high_positions[0]
        end = high_positions[0]

        for pos in high_positions[1:]:
            if pos == end + 1:
                end = pos
            else:
                if end - start + 1 >= min_region_length:
                    regions.append((start, end))    # save if long enough
                start = pos
                end = pos

        if end - start + 1 >= min_region_length:    # save the last region
            regions.append((start, end))

    return {
        "threshold": round(threshold, 5),
        "mean": round(mean, 5),
        "std": round(std, 5),
        "high_positions": high_positions,
        "cdr_regions": regions,
    }


def print_results(chain_name, results, entropy_values, chain_seqs):
    """
    Print a readable summary of what we found
    """

    regions = results["cdr_regions"]
    high_positions = results["high_positions"]

    print(f"\n Found {len(regions)} candidate CDR region(s):\n")

    if not regions:
        print(" None found. Try lowering STD_MULTIPLIER or MIN_REGION_LEN.")
        return

    seqs_list = list(chain_seqs.values())

    for i, (start, end) in enumerate(regions, 1):
        length = end - start + 1
        avg_entropy = sum(entropy_values[start:end+1]) / length

        print(f"  CDR candidate {i}: positions {start}-{end}  "
              f"(length={length}, avg_entropy={avg_entropy:.4f})")

        for header, seq in chain_seqs.items():
            parts = header.split("|")
            pdb_id = parts[0] if parts else header
            chain_desc = parts[2][:35] if len(parts) > 2 else ""
            region_seq = seq[start:end+1] if end+1 <= len(seq) else seq[start:]
            print(f"{pdb_id:12s} | {chain_desc:35s} | {region_seq}")
        print()

    # Show the 10 most variable positions so we can see what is driving the signal
    print(f" Top 10 most variable positions:")
    print(f"{'Pos':>5}  {'Entropy':>8}  AAs seen")
    print(f"{'-'*45}")
    top10 = sorted(enumerate(entropy_values), key=lambda x: -x[1])[:10]
    for pos, H in top10:
        aas = "".join(s[pos] for s in seqs_list if pos < len(s))
        flag = " <- high" if pos in high_positions else ""
        print(f"{pos:>5}  {H:>8.4f}  [{aas}]{flag}")



def save_results(chain_name, results, entropy_values, chain_seqs, output_path):
    """
    Save the full results to a text file. Includes all entropy values, regions found and sequences at each region
    """

    regions = results["cdr_regions"]
    high_positions = results["high_positions"]

    with open(output_path, "w") as f:
        f.write(f"CDR Detection Results - {chain_name} Chain\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sequences analysed: {len(chain_seqs)}\n")
        f.write(f"Positions: {len(entropy_values)}\n")
        f.write(f"Threshold: {results['threshold']}\n")
        f.write(f"Mean entropy: {results['mean']}\n")
        f.write(f"Std entropy: {results['std']}\n\n")

        f.write(f"Candidate CDR regions ({len(regions)} found):\n")
        f.write("-" * 60 + "\n")

        for i, (start, end) in enumerate(regions, 1):
            length = end - start + 1
            avg_H  = sum(entropy_values[start:end+1]) / length
            f.write(f"\nCDR candidate {i}: positions {start}-{end} "
                    f"(length={length}, avg_H={avg_H:.4f})\n")
            for header, seq in chain_seqs.items():
                parts = header.split("|")
                pdb_id = parts[0] if parts else header
                chain_desc = parts[2][:35] if len(parts) > 2 else ""
                region_seq = seq[start:end+1] if end+1 <= len(seq) else seq[start:]
                f.write(f" {pdb_id:12s} | {chain_desc:35s} | {region_seq}\n")

        f.write("\n\nFull entropy map:\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Pos':>5}  {'Entropy':>9}  Note\n")
        for pos, H in enumerate(entropy_values):
            flag = "  HIGH" if pos in high_positions else ""
            f.write(f"{pos:>5}  {H:>9.5f}{flag}\n")

    print(f"Saved: {output_path}")



def run_validation():
    """
    Test using the example from the הצעת מחקר:

        Pos:  0    1    2    3    4
        seq1: A    A    T    C    G
        seq2: A    T    T    C    G
        seq3: A    C    T    C    G
        seq4: A    G    T    G    G

    Position 0: all 'A' -> H should be 0.0  (no variability)
    Position 1: A, T, C, G -> H should be 2.0  (maximum for 4 options)
    Position 2: all 'T'  -> H should be 0.0  (no variability)
    """

    print("=" * 50)
    print("VALIDATION TEST")
    print("=" * 50)

    test_seqs = ["AATCG", "ATTCG", "ACTCG", "AGTGG"]
    entropy = calculate_entropy(test_seqs)

    print("Entropy per position:")
    for i, H in enumerate(entropy):
        print(f"  Position {i}: H = {H:.4f}")

    assert entropy[0] == 0.0,f"Expected H=0 at pos 0, got {entropy[0]}"
    assert abs(entropy[1]-2.0) < 0.001, f"Expected H=2.0 at pos 1, got {entropy[1]}"
    assert entropy[2] == 0.0, f"Expected H=0 at pos 2, got {entropy[2]}"
    assert entropy[1] == max(entropy), "Position 1 should have highest entropy"

    print("\n  All validation tests passed! v")
    print("=" * 50 + "\n")



def main():

    os.makedirs(RESULTS_DIR, exist_ok=True)

    #  make sure the math is correct on a known example
    run_validation()

    # Load all sequences from the FASTA file
    print(f"Loading sequences from: {FASTA_PATH}\n")
    all_seqs = parse_fasta(FASTA_PATH)

    if not all_seqs:
        print("No sequences found. Check the FASTA file path.")
        return

    print(f"\nTotal: {len(all_seqs)} sequences loaded.")

    # Split into heavy and light chains (different CDR structures)
    heavy = {h: s for h, s in all_seqs.items() if "heavy" in h.lower()}
    light = {h: s for h, s in all_seqs.items() if "light" in h.lower()}

    print(f"Heavy chains: {len(heavy)}")
    print(f"Light chains: {len(light)}")

    # Run analysis on each chain type
    for chain_name, chain_seqs in [("Heavy", heavy), ("Light", light)]:

        print(f"\n{'='*60}")
        print(f"  {chain_name.upper()} CHAIN ANALYSIS")
        print(f"{'='*60}")

        if len(chain_seqs) < 3:
            print(f"  Only {len(chain_seqs)} sequences — need at least 3. Skipping.")
            continue

        seqs_list = list(chain_seqs.values())

        # Calculate variability at each position
        entropy_values = calculate_entropy(seqs_list)

        # Find CDR candidate regions above the threshold
        results = find_cdr_regions(entropy_values, STD_MULTIPLIER, MIN_REGION_LEN)

        # Show results on screen
        print_results(chain_name, results, entropy_values, chain_seqs)

        # Save results to file
        out = os.path.join(RESULTS_DIR, f"{chain_name.lower()}_chain_results.txt")
        save_results(chain_name, results, entropy_values, chain_seqs, out)

        # Also compare with a stricter threshold so we can see the difference
        results_strict = find_cdr_regions(entropy_values, 2.0, MIN_REGION_LEN)
        print(f"\n  Threshold comparison:")
        print(f" mean + 1*std  ->  {len(results['cdr_regions'])} regions")
        print(f" mean + 2*std  ->  {len(results_strict['cdr_regions'])} regions")

    print(f"\nDone. Results saved to: {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
