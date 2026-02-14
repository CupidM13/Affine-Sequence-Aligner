import sys
import os
import datetime
import numpy as np

# ==============================================================================
# 1. FILE HANDLING
# ==============================================================================
def read_fasta(file_path):
    if not os.path.exists(file_path):
        print(f"\n[ERROR] File not found: '{file_path}'")
        sys.exit(1)

    sequence = []
    header_found = False
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if header_found: break 
                    header_found = True
                    continue
                sequence.append(line)
    except Exception as e:
        print(f"\n[ERROR] Could not read file: {e}")
        sys.exit(1)

    if not sequence:
        print(f"\n[ERROR] File '{file_path}' seems empty.")
        sys.exit(1)

    return "".join(sequence).upper()

def save_combined_report(seq1_path, full_seq1, seq2_path, full_seq2, 
                         res_global, res_local, params):
    """
    Saves a comprehensive report containing BOTH Global and Local alignments.
    """
    output_dir = "outputs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # --- SHORTENED FILENAME LOGIC ---
    # 1. Extract base name and take only first 10 chars
    name1 = os.path.splitext(os.path.basename(seq1_path))[0][:10]
    name2 = os.path.splitext(os.path.basename(seq2_path))[0][:10]
    
    # 2. Short timestamp: Day + Hour + Minute (e.g., 141230 for 14th 12:30)
    # This keeps it unique but much shorter than full YYYYMMDD...
    timestamp = datetime.datetime.now().strftime("%d%H%M")
    
    # Example result: geneA_geneB_141230.txt
    filename = f"{name1}_{name2}_{timestamp}.txt"
    filepath = os.path.join(output_dir, filename)

    # Unpack Results
    g_a1, g_a2, g_score, _ = res_global
    l_a1, l_a2, l_score, l_coords = res_local
    l_s1_start, l_s1_end, l_s2_start, l_s2_end = l_coords

    with open(filepath, "w") as f:
        # --- HEADER ---
        f.write("============================================================\n")
        f.write("          AFFINE GAP SEQUENCE ALIGNMENT REPORT\n")
        f.write("============================================================\n")
        f.write(f"Date:   {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"File 1: {os.path.basename(seq1_path)} (len: {len(full_seq1)})\n")
        f.write(f"File 2: {os.path.basename(seq2_path)} (len: {len(full_seq2)})\n")
        f.write("-" * 60 + "\n")
        f.write(f"PARAMETERS:\n")
        f.write(f"Match: {params['match']} | Mismatch: {params['mismatch']}\n")
        f.write(f"Gap Open (Sigma): {params['sigma']} | Gap Ext (Epsilon): {params['epsilon']}\n")
        f.write("=" * 60 + "\n\n")

        # --- SECTION 1: GLOBAL ALIGNMENT ---
        f.write("------------------------------------------------------------\n")
        f.write("  1. GLOBAL ALIGNMENT (Needleman-Wunsch / Gotoh)\n")
        f.write("------------------------------------------------------------\n")
        f.write(f"Strategy:    Aligns sequences from end-to-end.\n")
        f.write(f"Final Score: {g_score}\n\n")
        write_alignment_block(f, g_a1, g_a2)
        f.write("\n")

        # --- SECTION 2: LOCAL ALIGNMENT ---
        f.write("------------------------------------------------------------\n")
        f.write("  2. LOCAL ALIGNMENT (Smith-Waterman with Affine Gaps)\n")
        f.write("------------------------------------------------------------\n")
        f.write(f"Strategy:    Finds the best matching substring.\n")
        f.write(f"Final Score: {l_score}\n")
        f.write(f"Region 1:    Bases {l_s1_start+1} to {l_s1_end}\n")
        f.write(f"Region 2:    Bases {l_s2_start+1} to {l_s2_end}\n\n")
        write_alignment_block(f, l_a1, l_a2)
        f.write("\n")

        # --- SECTION 3: LOCAL CONTEXT VIEW ---
        f.write("------------------------------------------------------------\n")
        f.write("  3. LOCAL ALIGNMENT IN CONTEXT (Full Sequences)\n")
        f.write("------------------------------------------------------------\n")
        f.write("Showing the local match relative to the full sequence length.\n")
        f.write("Lowercase = Unaligned regions | Uppercase = Aligned Match\n\n")

        # Construct Context View
        prefix1 = full_seq1[:l_s1_start]
        suffix1 = full_seq1[l_s1_end:]
        prefix2 = full_seq2[:l_s2_start]
        suffix2 = full_seq2[l_s2_end:]

        # Pad prefixes to align matches vertically
        max_prefix = max(len(prefix1), len(prefix2))
        pad1 = " " * (max_prefix - len(prefix1))
        pad2 = " " * (max_prefix - len(prefix2))

        # Generate visual line
        core_vis = ""
        for a, b in zip(l_a1, l_a2):
            if a == b and a != '-': core_vis += "|"
            elif a == '-' or b == '-': core_vis += " "
            else: core_vis += "."
        
        full_vis = (" " * max_prefix) + core_vis
        
        # Build full strings
        global_s1 = pad1 + prefix1.lower() + l_a1 + suffix1.lower()
        global_s2 = pad2 + prefix2.lower() + l_a2 + suffix2.lower()

        write_alignment_block(f, global_s1, global_s2, full_vis)
        f.write("\n" + "="*60 + "\n")

    return filepath

def write_alignment_block(f, seq1, seq2, visual_line=None):
    """Helper to write wrapped alignment blocks"""
    block_width = 80
    generate_vis = (visual_line is None)

    for i in range(0, len(seq1), block_width):
        s1_slice = seq1[i:i+block_width]
        s2_slice = seq2[i:i+block_width]
        
        if generate_vis:
            vis_slice = ""
            for a, b in zip(s1_slice, s2_slice):
                if a == b and a != '-': vis_slice += "|"
                elif a == '-' or b == '-': vis_slice += " "
                else: vis_slice += "."
        else:
            if i < len(visual_line):
                end = min(i+block_width, len(visual_line))
                vis_slice = visual_line[i:end]
            else:
                vis_slice = ""

        f.write(f"Seq1: {s1_slice}\n")
        f.write(f"      {vis_slice}\n")
        f.write(f"Seq2: {s2_slice}\n\n")

# ==============================================================================
# 2. ALGORITHM (GOTOH / AFFINE GAP)
# ==============================================================================
def affine_gap_alignment(seq1, seq2, match, mismatch, sigma, epsilon, mode="global"):
    n = len(seq1)
    m = len(seq2)
    NEG_INF = -1e9
    
    # Matrices: Middle, Lower, Upper
    Middle = np.full((n + 1, m + 1), NEG_INF)
    Lower  = np.full((n + 1, m + 1), NEG_INF)
    Upper  = np.full((n + 1, m + 1), NEG_INF)

    Middle[0][0] = 0
    
    # Initialization
    if mode == "global":
        for i in range(1, n + 1):
            cost = sigma + (i - 1) * epsilon
            Lower[i][0], Middle[i][0] = cost, cost
        for j in range(1, m + 1):
            Upper[0][j], Middle[0][j] = cost, cost
    else:
        # Local: Borders 0
        for i in range(n + 1): Middle[i][0] = 0
        for j in range(m + 1): Middle[0][j] = 0

    # Fill Matrices
    max_score_local = -1
    max_pos_local = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_match = match if seq1[i-1] == seq2[j-1] else mismatch
            
            Lower[i][j] = max(Lower[i-1][j] + epsilon, Middle[i-1][j] + sigma)
            Upper[i][j] = max(Upper[i][j-1] + epsilon, Middle[i][j-1] + sigma)
            
            val = max(Middle[i-1][j-1] + score_match, Lower[i][j], Upper[i][j])
            
            if mode == "local":
                val = max(0, val)
                if val > max_score_local:
                    max_score_local = val
                    max_pos_local = (i, j)
            
            Middle[i][j] = val

    # Traceback
    align1, align2 = [], []
    
    if mode == "global":
        i, j = n, m
        scores = [Middle[n][m], Lower[n][m], Upper[n][m]]
        final_score = max(scores)
        state = scores.index(final_score)
        s1_end, s2_end = n, m
    else:
        i, j = max_pos_local
        final_score = max_score_local
        state = 0
        s1_end, s2_end = i, j

    while (mode == "global" and (i > 0 or j > 0)) or (mode == "local" and Middle[i][j] > 0):
        if state == 0: # Middle
            current = Middle[i][j]
            score_match = match if (i>0 and j>0 and seq1[i-1] == seq2[j-1]) else mismatch
            
            if i > 0 and j > 0 and current == Middle[i-1][j-1] + score_match:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1; j -= 1
            elif current == Lower[i][j]:
                state = 1
            elif current == Upper[i][j]:
                state = 2
            elif mode == "local" and current == 0:
                break
            else: 
                i -= 1; j -= 1
                 
        elif state == 1: # Lower
            align1.append(seq1[i-1])
            align2.append("-")
            if Lower[i][j] == Lower[i-1][j] + epsilon: state = 1
            else: state = 0
            i -= 1
            
        elif state == 2: # Upper
            align1.append("-")
            align2.append(seq2[j-1])
            if Upper[i][j] == Upper[i][j-1] + epsilon: state = 2
            else: state = 0
            j -= 1
            
    s1_start, s2_start = i, j
    return "".join(reversed(align1)), "".join(reversed(align2)), final_score, (s1_start, s1_end, s2_start, s2_end)

# ==============================================================================
# 3. MAIN (INPUT HANDLING)
# ==============================================================================
def get_user_input(prompt, default_val):
    user_val = input(f"{prompt} [Default: {default_val}]: ").strip()
    if not user_val: return default_val
    try: return float(user_val)
    except ValueError: return default_val

def main():
    print("\n" + "="*60)
    print("       AFFINE SEQUENCE ALIGNER (Global + Local)")
    print("="*60)
    
    if len(sys.argv) == 3:
        f1, f2 = sys.argv[1], sys.argv[2]
    else:
        f1 = input("File 1 (FASTA): ").strip()
        f2 = input("File 2 (FASTA): ").strip()

    s1 = read_fasta(f1)
    s2 = read_fasta(f2)
    print(f"Loaded: {len(s1)} bp vs {len(s2)} bp")

    # --- EXPLICIT WARNING FOR PENALTIES ---
    print("\n" + "-"*60)
    print("!!! IMPORTANT: SCORING PARAMETERS !!!")
    print("For Penalties (Mismatch, Gap Open, Gap Extend), you MUST enter NEGATIVE numbers.")
    print("Example: Enter '-5', NOT '5'. (Positive numbers are treated as rewards!)")
    print("-"*60)

    match = get_user_input("Match Score", 1.0)
    mismatch = get_user_input("Mismatch Score (Negative)", -1.0)
    sigma = get_user_input("Gap Open / Sigma (Negative)", -5.0)
    epsilon = get_user_input("Gap Extend / Epsilon (Negative)", -1.0)

    print("\nCalculating Global Alignment...")
    res_global = affine_gap_alignment(s1, s2, match, mismatch, sigma, epsilon, mode="global")
    
    print("Calculating Local Alignment...")
    res_local = affine_gap_alignment(s1, s2, match, mismatch, sigma, epsilon, mode="local")
    
    params = {'match': match, 'mismatch': mismatch, 'sigma': sigma, 'epsilon': epsilon}
    saved = save_combined_report(f1, s1, f2, s2, res_global, res_local, params)
    
    print("\n" + "="*60)
    print("SUCCESS!")
    print(f"Global Score: {res_global[2]}")
    print(f"Local Score:  {res_local[2]}")
    print(f"Report saved: {saved}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()