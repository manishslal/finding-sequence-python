# Finding_Sequence.py - Consensus Sequence Finder

This script analyzes a set of protein sequences from a FASTA file to find a consensus sequence and compare the original sequences against it.

## 1. Prerequisites

Before running the script, ensure you have the following installed:

* **Python 3:** The script is written in Python 3. Check your installation by running `python --version` or `python3 --version` in your terminal.

* **BioPython:** This library is used for parsing the FASTA file. Install it using pip:
    ```bash
    pip install biopython
    # or
    pip3 install biopython
    ```

* **(Optional) NCBI BLAST+:** If you want to use the (currently commented out) BLAST comparison feature, you need to have NCBI BLAST+ installed locally and configured in your system's PATH. You can download it from the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

## 2. Input File

* The script requires a single input file in **FASTA format**.
* This file should contain multiple protein sequences (amino acids).
* **Important:** The script assumes the sequences are **aligned and of the same length** for the consensus calculation (`finding_letters` function). If sequences have different lengths, the consensus calculation might be inaccurate, although the script includes a basic warning for this.
* The description line for each sequence in the FASTA file should ideally contain `Strain Name:` and `Protein Name:` tags for the script to extract the strain name correctly (though it has basic error handling if these are missing).

## 3. How to Run

1.  Open your terminal or command prompt.
2.  Navigate to the directory where you saved `Finding_Sequence.py` (e.g., `Finding-Sequence-Python`).
3.  Run the script using the following command structure:

    ```bash
    python Finding_Sequence.py <your_input_file.fasta>
    ```
    * Replace `<your_input_file.fasta>` with the actual path to your input FASTA file.

    **Example:**
    ```bash
    python Finding_Sequence.py my_protein_sequences.fasta
    ```

## 4. Output

The script will produce the following output:

* **Console Output:**
    * Confirmation of the input file being used.
    * The number of sequences parsed.
    * (Optional) Warnings if sequences have different lengths.
    * Information about the sequence that has the highest identity match to the generated consensus sequence (including percentage identity and match count).
    * A formatted table listing all input sequences (by Strain Name and Original ID part), sorted by their percentage identity match to the consensus sequence (highest first).
    * Confirmation that the consensus sequence was saved.
    * (If BLAST is enabled) Messages indicating the progress and completion of the `makeblastdb` and `blastp` commands.
    * The total time taken for the script to run.

* **Generated Files:**
    * `sample.fasta`: A new FASTA file containing the calculated consensus sequence. This file will be created in the same directory where you run the script.
    * (If BLAST is enabled) BLAST database files based on your input file name (e.g., `my_protein_sequences.fasta.phr`, `.pin`, `.psq`).
    * (If BLAST is enabled) `comparison_results.txt` (or the filename specified in `run_blast`): The output file containing the results of the BLASTp search comparing the consensus sequence against the database created from your input file.

## 5. Optional BLAST Comparison

The script includes a function (`run_blast`) to compare the generated consensus sequence (`sample.fasta`) against a BLAST database created from your original input file.

* **Enable:** To enable this, you need to uncomment the relevant lines near the end of the `if __name__ == "__main__":` block in `Finding_Sequence.py`:
    ```python
    # --- Optional: Run BLAST ---
    # Uncomment the following lines to run BLASTp
    # Requires BLAST+ installed and configured.
    blast_output_file = 'comparison_results.txt' # Choose your output filename
    run_blast(input_fasta_file, 'sample.fasta', blast_output_file)
    print(f"\nBLAST comparison attempted. Check '{blast_output_file}' for results (if successful).")
    ```
* **Requirements:** Ensure BLAST+ (`makeblastdb` and `blastp`) is installed and accessible from your command line.
* **Process:** When enabled, the script will first attempt to create a protein BLAST database from your input FASTA file. If successful, it will then run `blastp` using `sample.fasta` as the query against this newly created database. Results are saved to the specified output file (e.g., `comparison_results.txt`). Check the console output for any errors during the BLAST steps.

