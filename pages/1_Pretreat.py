import os
import io
import streamlit as st
import pandas as pd
from collections import Counter
from enum import Enum
import tempfile
from collections import defaultdict

# 导入要用于处理的函数
from utils import fasta_read, fasta_read2, fasta_read3, process_and_download, process_and_download2, get_column_names

def seq_length(file_path, min_length, max_length):
    """
    Filter sequences by length.

    Args:
        file_path (str): Path to the input FASTA file.
        min_length (int): Minimum sequence length.
        max_length (int): Maximum sequence length.

    Returns:
        str: Path to the output FASTA file containing filtered sequences.
    """
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        seq = seq.upper()
        if min_length < len(seq) < max_length:
            records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def name_change(file_path, info_path):
    """
    Replace sequence names according to the corresponding table.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to input Info file.

    Returns:
        str: Path to the output FASTA file with renamed sequences.
    """
    names, sequences = fasta_read2(file_path)
    df_info = pd.read_csv(info_path)
    old_name = df_info.iloc[:, 0].to_list()
    new_name = df_info.iloc[:, 1].to_list()
    name_dict = dict(zip(old_name, new_name))

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, sequence in zip(names, sequences):
                name_ = name_dict.get(name, name)
                f.write(f'>{name_}\n{sequence}\n')
    return out_path

class DeleteDuplicateFormat(str, Enum):
    DELETE = "delete"
    RENAME = "rename"

def handle_duplicate_ids(names, sequences, delete_format):
    """
    Handle duplicate IDs in sequences.

    Args:
        names (list): List of sequence IDs.
        sequences (list): List of sequences.
        delete_format (DeleteDuplicateFormat): Format for handling duplicates.

    Returns:
        tuple: Tuple containing new list of sequence IDs and sequences.
    """
    new_names = []
    new_sequences = []
    id_counts = {}
    for name, sequence in zip(names, sequences):
        if name in id_counts:
            if delete_format == DeleteDuplicateFormat.RENAME:
                id_counts[name] += 1
                new_name = f"{name}_{id_counts[name]}"
                new_names.append(new_name)
                new_sequences.append(sequence)
        else:
            id_counts[name] = 0
            new_names.append(name)
            new_sequences.append(sequence)
    return new_names, new_sequences

def id_duplicate(file_path, delete_format):
    """
    Handle duplicate sequence IDs in a FASTA file.

    Args:
        file_path (str): Path to the input FASTA file.
        delete_format (DeleteDuplicateFormat): Format for handling duplicates.

    Returns:
        str: Path to the output FASTA file with handled duplicate IDs.
    """
    names, sequences = fasta_read2(file_path)
    names, sequences = handle_duplicate_ids(names, sequences, delete_format)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, "w") as f:
            for name, sequence in zip(names, sequences):
                f.write(f">{name}\n{sequence}\n")
    return out_path

class ProcessFormat(str, Enum):
    DELETE = "Delete"
    REPLACE_N = "Replace to \"N\""
    REPLACE_HYPHEN = "Replace to \"-\""

def seq_standardize(file_path, standardize_format):
    """
    Process sequence characters to standardize them.

    Args:
        file_path (str): Path to the input FASTA file.
        standardize_format (ProcessFormat): Format for standardizing sequences.

    Returns:
        str: Path to the output FASTA file with standardized sequences.
    """
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        seq = seq.upper()
        if standardize_format == ProcessFormat.DELETE:
            seq = ''.join(c for c in seq if c in 'ATCG')
        elif standardize_format == ProcessFormat.REPLACE_N:
            seq = ''.join(c if c in 'ATCG' else 'N' for c in seq)
        elif standardize_format == ProcessFormat.REPLACE_HYPHEN:
            seq = ''.join(c if c in 'ATCG' else '-' for c in seq)

        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

class CaseFormat(str, Enum):
    UPPER = "upper"
    LOWER = "lower"

def seq_case(file_path, case_format):
    """
    Convert sequence case to upper or lower.

    Args:
        file_path (str): Path to the input FASTA file.
        case_format (CaseFormat): Format for case conversion.

    Returns:
        str: Path to the output FASTA file with converted case.
    """
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        if case_format == CaseFormat.UPPER:
            seq = seq.upper()
        elif case_format == CaseFormat.LOWER:
            seq = seq.lower()

        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def quality_control(file_path, qc_percentage):
    """
    Perform quality control on sequences based on ATCG content.

    Args:
        file_path (str): Path to the input FASTA file.
        qc_percentage (float): Minimum percentage of ATCG content.

    Returns:
        str: Path to the output FASTA file with sequences that meet the quality control criteria.
    """
    names, sequences = fasta_read2(file_path)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, sequence in zip(names, sequences):
                sequence = sequence.upper()
                seq_length = len(sequence)
                seq_count = Counter(sequence)
                atcg_count = seq_count["A"] + seq_count["T"] + seq_count["C"] + seq_count["G"]
                if atcg_count / seq_length >= qc_percentage:
                    f.write(f'>{name}\n{sequence}\n')
    return out_path

def reverse_sequence(file_path):
    """
    Reverse sequences in a FASTA file.

    Args:
        file_path (str): Path to the input FASTA file.

    Returns:
        str: Path to the output FASTA file with reversed sequences.
    """
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        seq = seq.upper()
        seq = seq[::-1]
        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def reverse_complement(file_path):
    """
    Reverse and complement sequences in a FASTA file.

    Args:
        file_path (str): Path to the input FASTA file.

    Returns:
        str: Path to the output FASTA file with reversed and complemented sequences.
    """
    sequences = fasta_read(file_path)
    records = []
    for name, seq in sequences.items():
        seq = seq.upper()
        seq = seq[::-1]
        seq = seq.translate(str.maketrans("ATCG", "TAGC"))
        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def merge_multi_fasta(fasta_files):
    """
    Merge multiple FASTA files into a single FASTA file.

    Args:
        fasta_files (List[str]): List of paths to FASTA files.

    Returns:
        str: Path to the output FASTA file containing merged sequences.
    """
    records = []
    for fasta_file in fasta_files:
        sequences = fasta_read3(fasta_file)
        for name, seq in sequences.items():
            records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def main():
    """
    Main function to run the Streamlit app for sequence pretreatment.
    """
    st.title("Sequence Pretreatment")

    selected_option = st.radio("Select an option", ["Filtered by length", 
                                                    "Replace sequence id", 
                                                    "Handle duplicate IDs", 
                                                    "Process sequence characters", 
                                                    "Convert sequence case", 
                                                    "Quality control",
                                                    "Reverse sequence",
                                                    "Reverse and complement sequence",
                                                    "Merge multiple FASTA files"])

    if selected_option == "Filtered by length":
        min_length = st.number_input("Min length", step=1)
        max_length = st.number_input("Max length", step=1)
        file_path = st.file_uploader("Upload a Fasta file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: seq_length(x, min_length, max_length), "length_filtered.fasta")

    elif selected_option == "Replace sequence id":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        info_path = st.file_uploader("Upload an Info file", type=["csv", "table", "txt", "xlsx"])
        if file_path and info_path:
            process_and_download2(file_path, info_path, lambda x, y: name_change(x, y), "renamed_sequences.fasta")

    elif selected_option == "Handle duplicate IDs":
        delete_format = st.selectbox("Delete Duplicate Format", ["rename", "delete"])
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: id_duplicate(x, DeleteDuplicateFormat(delete_format)), "handled_duplicate_ids.fasta")

    elif selected_option == "Process sequence characters":
        standardize_format = st.selectbox("Standardize Format", ["Delete", "Replace to \"N\"", "Replace to \"-\""])
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: seq_standardize(x, ProcessFormat(standardize_format)), "standardized_sequences.fasta")

    elif selected_option == "Convert sequence case":
        case_format = st.selectbox("Case Format", ["upper", "lower"])
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: seq_case(x, CaseFormat(case_format)), "converted_case_sequences.fasta")

    elif selected_option == "Quality control":
        qc_percentage = st.slider("QC Percentage", min_value=0.0, max_value=1.0, step=0.01)
        file_path = st.file_uploader("Upload a Fasta file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: quality_control(x, qc_percentage), "quality_controlled_sequences.fasta")
    
    elif selected_option == "Reverse sequence":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: reverse_sequence(x), "reversed_sequences.fasta")
    
    elif selected_option == "Reverse and complement sequence":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        if file_path:
            process_and_download(file_path, lambda x: reverse_complement(x), "reversed_and_complemented_sequences.fasta")

    elif selected_option == "Merge multiple FASTA files":
        fasta_files = st.file_uploader("Upload multiple FASTA files", type=["fasta", "fas", "fa"], accept_multiple_files=True)
        if fasta_files:
            if st.button("Run"):
                out_file_path = merge_multi_fasta(fasta_files)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name="merged.fasta",
                    mime="text/x-fasta"
                )
                os.remove(out_file_path)
    


if __name__ == "__main__":
    main()
