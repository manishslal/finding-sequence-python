import collections
from Bio import SeqIO
import sys
import time
import math
import subprocess as sp
import os # Import os module for file existence check

# Initialize lists
list1, list2, list3, list3a, final_list, new_list = [], [], [], [], [], []
# list1: Strain names + years
# list2: Sequences
# list3: Strings of amino acids at each position
# list3a: Original sequence IDs (up to '|')
# final_list: Consensus sequence (list of chars)
# new_list: Sequence identity percentages

# --- Function Definitions ---

def getting_filename():
    """
    Gets the FASTA filename from command-line arguments.
    Exits if no filename is provided.
    """
    if len(sys.argv) > 1:
        filename_path = sys.argv[1]
        # Check if the file exists before trying to open
        if os.path.exists(filename_path):
            print(f"Using input file: {filename_path}")
            return filename_path
        else:
            print(f"Error: File not found at '{filename_path}'")
            sys.exit(1) # Exit if file doesn't exist
    else:
        print('Error: No input FASTA file name provided.')
        print('Usage: python Finding_Sequence.py <your_fasta_file.fasta>')
        sys.exit(1) # Exit if no filename is given

def getting_aa(filename_path, list_one, list_two, list_three_a):
    """
    Reads through the FASTA file using Biopython, extracts sequence IDs,
    derived strain names, and sequences into respective lists.
    Args:
        filename_path (str): Path to the input FASTA file.
        list_one (list): List to store derived strain names.
        list_two (list): List to store sequences.
        list_three_a (list): List to store original sequence IDs (partial).
    """
    i = 0
    try:
        # Use SeqIO.parse which returns an iterator
        for record in SeqIO.parse(filename_path, "fasta"):
            ids = record.id   # Takes the ID of the sequence
            seq = record.seq  # Takes the sequence object
            description = record.description

            # Taking in just the name of the Strain and its isolation year for future reference
            # Handle cases where 'Strain Name:' might be missing
            name = "Unknown Strain" # Default name
            try:
                start = description.find('Strain Name:') + 12
                end = description.find('Protein Name:')
                if start != 11 and end != -1: # Check if find() was successful
                     name = description[start:end].strip()
                elif start != 11: # If only start is found
                     name = description[start:].strip()
                # Add more robust parsing if needed based on actual description formats
            except Exception as e:
                print(f"Warning: Could not parse Strain Name from description: {description} - Error: {e}")


            list_one.append(name)
            list_two.append(str(seq)) # Store sequence as string
            # Store original ID part safely
            id_part = ids
            if '|' in ids:
                id_part = ids[:ids.index('|')]
            list_three_a.append(id_part)
            i += 1
        print(f"Successfully parsed {i} sequences.")
        # Check if sequences were actually read
        if not list_two:
             print("Error: No sequences found in the input file.")
             sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Input file '{filename_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during FASTA parsing: {e}")
        sys.exit(1)


def finding_letters(list_of_seqs, list_letters_at_pos):
    """
    Extracts individual letters at each position across all sequences
    and creates strings of these letters. Assumes sequences are aligned
    and have the same length.
    Args:
        list_of_seqs (list): List containing sequence strings.
        list_letters_at_pos (list): List to store strings of letters at each position.
    """
    if not list_of_seqs:
        print("Error: No sequences available to find letters.")
        return

    # Check if all sequences have the same length (crucial for this method)
    first_len = len(list_of_seqs[0])
    if not all(len(seq) == first_len for seq in list_of_seqs):
        print("Warning: Input sequences have different lengths. Consensus calculation might be inaccurate.")
        # Decide how to handle: truncate, pad, or error out?
        # For now, proceed using the length of the first sequence. Adjust if needed.

    j = 0
    while j < first_len:    # Loop through each position
        k = 0
        letter_string = ''
        while k < len(list_of_seqs):   # Loop through each sequence
            # Add boundary check in case of differing lengths despite warning
            if j < len(list_of_seqs[k]):
                letter_string = list_of_seqs[k][j] + letter_string # Build string of letters at position j
            k += 1
        j += 1
        list_letters_at_pos.append(letter_string) # Add the string for this position

def finding_frequency(list_letters_at_pos, consensus_list):
    """
    Finds the most common letter (amino acid/nucleotide) at each position.
    Args:
        list_letters_at_pos (list): List of strings, each containing letters at one position.
        consensus_list (list): List to store the most common letter for each position.
    """
    m = 0
    while m < len(list_letters_at_pos):   # loop that reads through the list of position strings
        string_at_pos = list_letters_at_pos[m]
        if not string_at_pos: # Handle empty string case
            consensus_list.append('-') # Or 'N', or 'X' depending on context
            print(f"Warning: No characters found at position {m}")
        else:
            # Find the most frequent letter using Counter
            most_common_tuple = collections.Counter(string_at_pos).most_common(1)[0]
            consensus_list.append(most_common_tuple[0]) # Add the most common letter
        m += 1

def finding_sequence(list_of_sequences, consensus_sequence_list, list_strain_names, list_identity_percentages, list_original_ids):
    """
    Compares all original sequences to the generated consensus sequence, calculates
    percentage identity, and prints results.
    Args:
        list_of_sequences (list): List of original sequence strings.
        consensus_sequence_list (list): Consensus sequence as a list of characters.
        list_strain_names (list): List of derived strain names.
        list_identity_percentages (list): List to store calculated identity percentages.
        list_original_ids (list): List of original sequence IDs (partial).
    """
    if not list_of_sequences or not consensus_sequence_list:
        print("Error: Cannot compare sequences. Input lists are empty.")
        return

    consensus_len = len(consensus_sequence_list)
    consensus_sequence_str = "".join(consensus_sequence_list) # Join for easier comparison if needed

    print("\nComparing original sequences to consensus sequence...")

    highest_identity = -1.0
    best_match_index = -1
    total_compared_overall = 0
    total_matches_overall = 0 # For overall average if needed

    for idx, seq_str in enumerate(list_of_sequences):
        match_count = 0
        # Compare only up to the length of the shorter sequence (consensus or current seq)
        comparison_len = min(len(seq_str), consensus_len)
        if comparison_len == 0:
             list_identity_percentages.append(0.0)
             continue # Skip empty sequences

        for i in range(comparison_len):
            if seq_str[i] == consensus_sequence_list[i]:
                match_count += 1

        # Calculate percentage based on the comparison length
        identity_percentage = (float(match_count) / float(comparison_len)) * 100.0
        list_identity_percentages.append(identity_percentage)

        # Track overall stats if needed
        total_matches_overall += match_count
        total_compared_overall += comparison_len

        # Find the sequence with the highest identity
        if identity_percentage > highest_identity:
            highest_identity = identity_percentage
            best_match_index = idx

    # Print the best match result
    if best_match_index != -1:
        # Calculate matches/total for the specific best match sequence
        best_match_seq = list_of_sequences[best_match_index]
        best_match_len = min(len(best_match_seq), consensus_len)
        best_match_count = 0
        for i in range(best_match_len):
             if best_match_seq[i] == consensus_sequence_list[i]:
                 best_match_count += 1

        print(f"\nHighest sequence identity match found:")
        print(f"  Sequence: {list_strain_names[best_match_index]} ({list_original_ids[best_match_index]})")
        print(f"  Identity: {highest_identity:.2f}% ({best_match_count}/{best_match_len} matches)")
    else:
        print("\nCould not determine the best matching sequence.")


    # Print the sorted table of all sequences and their identity matches
    print("\n--- Sequence Identity Comparison ---")
    print(f"{'Strain Name':<30} | {'Original ID':<15} | {'Identity Match (%)':>18}") # Adjusted column widths
    print("-" * 70) # Separator line

    # Create pairs of (identity, index) to sort, handling potential length mismatch
    identity_pairs = []
    if len(list_identity_percentages) == len(list_strain_names) == len(list_original_ids):
         identity_pairs = sorted([(list_identity_percentages[i], i) for i in range(len(list_identity_percentages))], reverse=True)
    else:
         print("Warning: List length mismatch, cannot print detailed comparison table.")
         return # Exit if lists don't match

    for identity, index in identity_pairs:
        strain_name = list_strain_names[index]
        original_id = list_original_ids[index]
        print(f"{strain_name:<30} | {original_id:<15} | {identity:>18.2f}")


def save_sequence(consensus_sequence_list):
    """
    Saves the generated consensus sequence into a FASTA file named 'sample.fasta'.
    Args:
        consensus_sequence_list (list): Consensus sequence as a list of characters.
    """
    output_filename = 'sample.fasta'
    print(f"\nSaving consensus sequence to {output_filename}...")
    try:
        with open(output_filename, 'w') as output:
            # Basic FASTA header (customize as needed)
            output.write(">Consensus_Sequence|Generated_by_Finding_Sequence.py\n")

            consensus_string = "".join(consensus_sequence_list)
            # Write sequence with line breaks (e.g., every 70 characters)
            line_length = 70
            for i in range(0, len(consensus_string), line_length):
                output.write(consensus_string[i:i+line_length] + '\n')
            print(f"Consensus sequence saved successfully.")
            # Optionally print the sequence to console as well
            # print('\nConsensus Sequence:')
            # print(consensus_string)

    except IOError as e:
        print(f"Error: Could not write to file {output_filename}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred during sequence saving: {e}")


def run_blast(subject_db_path, query_fasta_path, output_filename):
    """
    Runs a local BLASTp search using the provided paths.
    Requires BLAST+ to be installed and configured in the system PATH or provide full path.
    Args:
        subject_db_path (str): Path to the BLAST database (created from the input fasta).
        query_fasta_path (str): Path to the query FASTA file (e.g., 'sample.fasta').
        output_filename (str): Name for the BLAST output file.
    """
    # --- IMPORTANT: Create BLAST DB First! ---
    # BLASTp needs a database, not just the FASTA file directly as subject.
    # You need to run makeblastdb first.
    db_creation_cmd = [
        'makeblastdb',
        '-in', subject_db_path, # Input FASTA used to create DB
        '-dbtype', 'prot',      # Database type is protein
        '-out', subject_db_path # Base name for database files
    ]
    print(f"\nAttempting to create BLAST database from {subject_db_path}...")
    try:
        # Run makeblastdb and wait for it to complete
        process_db = sp.run(db_creation_cmd, check=True, capture_output=True, text=True)
        print("BLAST database created successfully.")
        # print(process_db.stdout) # Optional: print stdout
        # print(process_db.stderr) # Optional: print stderr

        # --- Now Run BLASTp ---
        blast_cmd = [
            # Provide full path if needed, e.g., '/usr/local/ncbi/blast/bin/blastp'
            'blastp',
            '-db', subject_db_path,    # Use the created database base name
            '-query', query_fasta_path, # The consensus sequence file
            '-out', output_filename,    # Where to save results
            '-outfmt', '6'             # Example: Tabular output format (customize as needed)
            # Add other BLAST options here if desired (e.g., -evalue, -max_target_seqs)
        ]
        print(f"Running BLASTp: {' '.join(blast_cmd)}")
        # Use Popen if you want it to run in the background, or run() to wait for completion
        # Using run() is often better for scripts unless background processing is needed.
        process_blast = sp.run(blast_cmd, check=True, capture_output=True, text=True)
        print(f"BLASTp search completed. Results saved to {output_filename}")
        # print(process_blast.stdout) # Optional: print stdout from blastp
        # print(process_blast.stderr) # Optional: print stderr from blastp

    except FileNotFoundError as e:
         print(f"Error: BLAST command not found. Is BLAST+ installed and in your system PATH? ({e})")
    except sp.CalledProcessError as e:
         print(f"Error during BLAST execution (return code {e.returncode}):")
         print(f"  Command: {' '.join(e.cmd)}")
         print(f"  Stderr: {e.stderr}")
    except Exception as e:
         print(f"An unexpected error occurred during BLAST execution: {e}")


# --- Main Execution ---
if __name__ == "__main__":
    start_time = time.time()

    # Get filename from command line
    input_fasta_file = getting_filename()

    # Run the analysis steps
    getting_aa(input_fasta_file, list1, list2, list3a)
    finding_letters(list2, list3)
    finding_frequency(list3, final_list)

    # Check if consensus sequence was generated
    if final_list:
        # The final_list is converted into a string, so that it can be printed/saved
        produced_sequence_list = final_list # Keep as list for comparison function
        save_sequence(produced_sequence_list)
        time.sleep(0.3) # Keep small delay if desired
        finding_sequence(list2, produced_sequence_list, list1, new_list, list3a)

        # --- Optional: Run BLAST ---
        # Uncomment the following lines to run BLASTp
        # Requires BLAST+ installed and configured.
        # blast_output_file = 'comparison_results.txt'
        # run_blast(input_fasta_file, 'sample.fasta', blast_output_file)
        # print(f"\nBLAST comparison attempted. Check '{blast_output_file}' for results (if successful).")

    else:
        print("\nConsensus sequence could not be generated. Skipping comparison and saving.")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"\nScript finished in {elapsed_time:.2f} seconds.")
