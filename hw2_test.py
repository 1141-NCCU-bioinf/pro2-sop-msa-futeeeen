import pandas as pd

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of sequences.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        list: A list of strings, where each string is a sequence.
    """
    sequences = []
    current_sequence = ""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_sequence:
                        sequences.append(current_sequence)
                    current_sequence = ""
                else:
                    current_sequence += line
            if current_sequence:
                sequences.append(current_sequence)
    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return []
    return sequences

def read_score_matrix(file_path):
    """
    Reads a scoring matrix file (like BLOSUM or PAM).

    Args:
        file_path (str): The path to the scoring matrix file.

    Returns:
        pandas.DataFrame: A DataFrame representing the scoring matrix.
    """
    try:
        # Reads a whitespace-delimited file, ignoring lines starting with '#'
        # Using sep='\s+' instead of the deprecated delim_whitespace=True
        matrix = pd.read_csv(file_path, sep=r'\s+', comment='#')
        return matrix
    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return None

def calculate_SoP(input_path, score_path, gopen, gextend):
    """
    Calculates the Sum-of-Pairs (SoP) score for a multiple sequence alignment.

    Args:
        input_path (str): Path to the input FASTA file with the alignment.
        
        score_path (str): Path to the scoring matrix file.
        gopen (int): The penalty for opening a gap.
        gextend (int): The penalty for extending a gap.

    Returns:
        int or None: The total Sum-of-Pairs score, or None if an error occurs.
    """
    sequences = read_fasta(input_path)
    score_matrix = read_score_matrix(score_path)

    if not sequences or score_matrix is None:
        print("Could not read input files. Aborting.")
        return None

    num_sequences = len(sequences)
    if num_sequences < 2:
        return 0  # Score is 0 if there are fewer than 2 sequences.

    total_sop_score = 0

    # Iterate over all unique pairs of sequences
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            seq1 = sequences[i]
            seq2 = sequences[j]
            
            pair_score = 0
            is_gap_in_seq1 = False
            is_gap_in_seq2 = False

            # Calculate the score for the current pair
            for char1, char2 in zip(seq1, seq2):
                each_score = 0
                if char1 != '-' and char2 != '-':
                    # Case 1: Both are amino acids (match or mismatch)
                    try:
                        each_score += score_matrix.loc[char1, char2]
                    except KeyError:
                        print(f"Warning: Character pair ('{char1}', '{char2}') not found in score matrix. Treating as 0.")
                    is_gap_in_seq1 = False
                    is_gap_in_seq2 = False
                elif char1 == '-' and char2 != '-':
                    # Case 2: Gap in sequence 1
                    each_score += gextend if is_gap_in_seq1 else gopen
                    is_gap_in_seq1 = True
                    is_gap_in_seq2 = False
                elif char1 != '-' and char2 == '-':
                    # Case 3: Gap in sequence 2
                    each_score += gextend if is_gap_in_seq2 else gopen
                    is_gap_in_seq2 = True
                    is_gap_in_seq1 = False
                else:  # Both are gaps ('-')
                    # Case 4: Gap-gap alignment has no cost
                    each_score += gextend if (is_gap_in_seq2 or is_gap_in_seq1) else gopen
                    is_gap_in_seq1 = True
                    is_gap_in_seq2 = True
                print(f"[{char1}],[{char2}] score : {each_score}")
                pair_score += each_score
            print(f"score : {pair_score}")
            total_sop_score += pair_score

    return total_sop_score

if __name__ == '__main__':
    # --- Parameters ---
    # Make sure these files are in the same directory as the script,
    # or provide the full path.
    fasta_file = 'examples/test1.fasta'
    score_file = 'examples/pam250.txt'
    
    # Gap penalties should typically be negative
    gap_open_penalty = -10
    gap_extend_penalty = -2

    # --- Calculation ---
    print(f"Calculating SoP score for '{fasta_file}' using '{score_file}'...")
    print(f"Gap Open Penalty: {gap_open_penalty}, Gap Extend Penalty: {gap_extend_penalty}\n")
    
    sop_score = calculate_SoP(fasta_file, score_file, gap_open_penalty, gap_extend_penalty)

    # --- Output ---
    if sop_score is not None:
        print(f"The final Sum-of-Pairs score is: {sop_score}")
