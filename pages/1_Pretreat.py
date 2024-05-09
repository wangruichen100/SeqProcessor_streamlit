import os
import io
import streamlit as st
import pandas as pd
from Bio import SeqIO
from collections import Counter
from enum import Enum
import tempfile

# 导入要用于处理的函数
from utils import fasta_read, fasta_read2


def name_change(file_path, info_path):
    """
    Replace sequence names according to the corresponding table.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to input Info file.
        out_path (str): Path to the output file.
        info_format (info_format): Format for Info file format.

    Returns:
        None
    """
    names, sequences = fasta_read2(file_path)
    
    df_info = pd.read_csv(info_path)
    
    old_name = df_info.iloc[:, 0].to_list()
    new_name = df_info.iloc[:, 1].to_list()
    name_dict = dict(zip(old_name,new_name))
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, sequence in zip(names, sequences):
                name_ = name_dict.get(name, name)
                f.write(f'>{name_}\n{sequence}\n')

    return out_path

# 处理重复的 ID
class DeleteDuplicateFormat(str, Enum):
    f1 = "delete"  # "Delete duplicate sequences."
    f2 = "rename"  # "Rename duplicate sequence names."
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
            if delete_format == DeleteDuplicateFormat.f2:
                # If duplicate IDs should be deleted, add a suffix to the name
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
    names, sequences = fasta_read2(file_path)
    names, sequences = handle_duplicate_ids(names, sequences, delete_format)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, "w") as f:
            for name, sequence in zip(names, sequences):
                f.write(f">{name}\n{sequence}\n")
    return out_path

# 处理序列中的字符
class ProcessFormat(str, Enum):
    f1 = "Delete"  # "Delete characters outside of ATCG."
    f2 = "Replace to \"N\""  # "Replace characters outside of ATCG with N."
    f3 = "Replace to \"-\""  # "Replace characters outside of ATCG with -."
def seq_standardize(file_path, standardize_format):
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        seq = seq.upper()
        if standardize_format == ProcessFormat.f1:
            seq = ''.join(c for c in seq if c in 'ATCGatcg')
        elif standardize_format == ProcessFormat.f2:
            seq = ''.join(c if c in 'ATCGatcg' else 'N' for c in seq)
        elif standardize_format == ProcessFormat.f3:
            seq = ''.join(c if c in 'ATCGatcg' else '-' for c in seq)

        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

# 转换序列的大小写
class CaseFormat(str, Enum):
    f1 = "upper"  # "Convert to uppercase."
    f2 = "lower"  # "Convert to lowercase."
def seq_case(file_path, case_format):
    sequences = fasta_read(file_path)
    records = []

    for name, seq in sequences.items():
        if case_format == CaseFormat.f1:
            seq = seq.upper()
        elif case_format == CaseFormat.f2:
            seq = seq.lower()

        records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

# 过滤低质量序列
def quality_control(file_path, qc_percentage):
    names, sequences = fasta_read2(file_path)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, sequence in zip(names, sequences):
                sequence = sequence.upper()
                seq_length = len(sequence)
                seq_count = Counter(sequence)
                atcg_count = seq_count["A"]+seq_count["T"]+seq_count["C"]+seq_count["G"]
                if atcg_count/seq_length >= qc_percentage:
                    f.write(f'>{name}\n{sequence}\n')
    return out_path

# 创建 Streamlit 应用程序
def process_and_download(file_path, process_function, output_filename):
    if file_path:
        if st.button("Process"):
            # Save uploaded file to a temporary location
            temp_file_path = tempfile.NamedTemporaryFile(delete=True, suffix=".fasta").name
            with open(temp_file_path, "wb") as f:
                f.write(file_path.getvalue())
            
            out_file_path = process_function(temp_file_path)
            st.success("File processed successfully!")
            with open(out_file_path, "r") as f:
                processed_content = f.read()
            st.download_button(
                label="Download Processed File",
                data=processed_content,
                file_name=output_filename,  # Use the defined output filename here
                mime="application/octet-stream"
            )
            os.remove(out_file_path)
            os.remove(temp_file_path)

def process_and_download2(file_path, info_path, output_filename):
    if file_path and info_path:
        if st.button("Process"):
            # Save uploaded file to a temporary location
            temp_fas_file_path = tempfile.NamedTemporaryFile(delete=True, suffix=".fasta").name
            with open(temp_fas_file_path, "wb") as f:
                f.write(file_path.getvalue())

            # Determine the format of the info file based on its extension
            info_format = info_path.name.split(".")[-1]
            temp_info_file_path = tempfile.NamedTemporaryFile(delete=True).name
            if info_format == "csv":
                df = pd.read_csv(info_path)
                df.to_csv(temp_info_file_path, index=False)
            elif info_format == "txt" or info_format == "tsv" or info_format == "tab":
                df = pd.read_csv(info_path, sep="\t")
                df.to_csv(temp_info_file_path, index=False)
            elif info_format == "xlsx":
                df = pd.read_excel(info_path, engine='openpyxl')
                df.to_csv(temp_info_file_path, index=False)

            out_file_path = name_change(temp_fas_file_path, temp_info_file_path)
            st.success("File processed successfully!")
            with open(out_file_path, "r") as f:
                processed_content = f.read()
            st.download_button(
                label="Download Processed File",
                data=processed_content,
                file_name=output_filename,  # Use the defined output filename here
                mime="application/octet-stream"
            )


def main():
    st.title("Sequence Pretreatment")

    # 创建顶部选项
    selected_option = st.radio("Select an option", ["Replace sequence names", "Handle Duplicate IDs", "Process Sequence Characters", "Convert Sequence Case", "Quality Control"])

    if selected_option == "Replace sequence names":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas"])
        info_path = st.file_uploader("Upload an Info file", type=["csv", "table", "xlsx"])
        process_and_download2(file_path, info_path, "renamed_sequences.fasta")  # Define the output filename here

    elif selected_option == "Handle Duplicate IDs":
        delete_format = st.selectbox("Delete Duplicate Format", ["rename","delete"])

        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas"])
        process_and_download(file_path, lambda x: id_duplicate(x, delete_format), "handled_duplicate_ids.fasta")  # Define the output filename here

    elif selected_option == "Process Sequence Characters":
        standardize_format = st.selectbox("Standardize Format", ["Delete", "Replace to \"N\"", "Replace to \"-\""])

        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas"])
        process_and_download(file_path, lambda x: seq_standardize(x, standardize_format), "standardized_sequences.fasta")  # Define the output filename here

    elif selected_option == "Convert Sequence Case":
        case_format = st.selectbox("Case Format", ["upper", "lower"])

        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas"])
        process_and_download(file_path, lambda x: seq_case(x, case_format), "converted_case_sequences.fasta")  # Define the output filename here

    elif selected_option == "Quality Control":
        qc_percentage = st.slider("QC Percentage", min_value=0.0, max_value=1.0, step=0.01)
        file_path = st.file_uploader("Upload a Fasta file", type=["fasta", "fas"])
        process_and_download(file_path, lambda x: quality_control(x, qc_percentage), "quality_controlled_sequences.fasta")  # Define the output filename here

if __name__ == "__main__":
    main()
