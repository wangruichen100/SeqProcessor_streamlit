import os
import io
import streamlit as st
import pandas as pd
import tempfile
import random
import shutil

# 导入要用于处理的函数
from utils import fasta_read, fasta_read2, get_column_names

def extract_sequence(file_path, info_path):
    """
    Extract sequences from a FASTA file based on IDs in an info file.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to the info file containing sequence IDs.

    Returns:
        str: Path to the output FASTA file containing extracted sequences.
    """
    sequences = fasta_read(file_path)
    df_info = pd.read_csv(info_path)
    require_id = df_info.iloc[:, 0].to_list()

    records = []

    for name, seq in sequences.items():
        seq = seq.upper()
        for id_ in require_id:
            if id_ in name:
                records.append({"Name": name, "Sequence": seq})

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for record in records:
                f.write(f'>{record["Name"]}\n{record["Sequence"]}\n')
    return out_path

def random_rate_sequence(file_path, random_rate):
    """
    Extract a random subset of sequences from a FASTA file based on a specified rate.

    Args:
        file_path (str): Path to the input FASTA file.
        random_rate (float): The rate of sequences to be randomly selected.

    Returns:
        str: Path to the output FASTA file containing the random subset of sequences.
    """
    sequences = fasta_read(file_path)
    num_items_to_select = int(len(sequences) * random_rate)
    sequence_list = list(sequences.items())
    sequences_random = random.sample(sequence_list, num_items_to_select)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, seq in sequences_random:
                f.write(f'>{name}\n{seq}\n')
    return out_path

def random_number_sequence(file_path, num_sequences):
    """
    Extract a specified number of random sequences from a FASTA file.

    Args:
        file_path (str): Path to the input FASTA file.
        num_sequences (int): The number of sequences to be randomly selected.

    Returns:
        str: Path to the output FASTA file containing the random subset of sequences.
    """
    sequences = fasta_read(file_path)
    sequence_names = list(sequences.keys())
    selected_names = random.sample(sequence_names, num_sequences)
    selected_sequences = {name: sequences[name] for name in selected_names}

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, seq in selected_sequences.items():
                f.write(f'>{name}\n{seq}\n')
    return out_path

def group_sequence(file_path, info_path, id_column, type_column):
    """
    Extract sequences from a FASTA file and group them based on information in an info file.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to the info file containing grouping information.
        id_column (str): Column name for sequence IDs in the info file.
        type_column (str): Column name for grouping criteria in the info file.

    Returns:
        str: Path to the output ZIP file containing grouped sequences.
    """
    seq_record = fasta_read(file_path)
    df = pd.read_excel(info_path)
    seq_types = sorted(set(df[type_column].to_list()))
    output_folder = tempfile.mkdtemp()

    for seq_type in seq_types:
        df_type = df[df[type_column] == seq_type]
        with open(os.path.join(output_folder, f"{seq_type}.fas"), "w") as f:
            for name in df_type[id_column].to_list():
                f.write(f">{name}\n{seq_record[name]}\n")

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename


def extract_gene(file_path, info_path, gene_column, start_column, end_column):
    names, sequences = fasta_read2(file_path)
    # Read the gene information from the Excel file
    df = pd.read_excel(info_path)
    genes = df[gene_column].to_list()
    starts = df[start_column].to_list()
    ends = df[end_column].to_list()

    output_folder = tempfile.mkdtemp()

    # Iterate over each gene, start, and end position
    for gene, start, end in zip(genes, starts, ends):
        print(gene)

        if start > end:
            start, end = end, start
        # Open the output file for the current gene
        batch_filename = os.path.join(output_folder, f"{gene}.fasta")
        with open(batch_filename, "w") as f:
            # Read sequences from the FASTA file
            for name, seq in zip(names, sequences):
                # Extract the subsequence from start to end (1-based indexing)
                subsequence = seq[start-1:end]
                # Write the subsequence to the output file
                f.write(f">{name}\n{subsequence}\n")

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename

def split_fasta(file_path, sequences_per_file):
    """
    Split a FASTA file into multiple files, each containing a specified number of sequences.

    Args:
        file_path (str): Path to the input FASTA file.
        sequences_per_file (int): Number of sequences per output file.

    Returns:
        str: Path to the output ZIP file containing the split FASTA files.
    """
    sequences = fasta_read(file_path)
    sequence_items = list(sequences.items())
    output_folder = tempfile.mkdtemp()

    for i in range(0, len(sequence_items), sequences_per_file):
        batch = sequence_items[i:i + sequences_per_file]
        batch_filename = os.path.join(output_folder, f"batch_{i // sequences_per_file + 1}.fasta")
        with open(batch_filename, 'w') as f:
            for name, seq in batch:
                f.write(f'>{name}\n{seq}\n')

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename


def main():
    """
    Main function to run the Streamlit app for sequence extraction.
    """
    st.title("Sequence Extraction")

    selected_option = st.radio("Select an option", ["Extract sequence by id",
                                                    "Extract sequence by random rate",
                                                    "Extract sequence by random number",
                                                    "Extract sequence by group",
                                                    "Extract sequence by gene site",
                                                    "Split FASTA file into multiple files"])
    
    if selected_option == "Extract sequence by id":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        info_path = st.file_uploader("Upload an Info file", type=["csv", "table", "txt", "xlsx"])
        
        if file_path and info_path:
            if st.button("Run"):
                temp_file_path = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                out_file_path = extract_sequence(temp_file_path, info_path)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name="require_seq.fasta",
                    mime="application/octet-stream"
                )
                os.remove(out_file_path)
                os.remove(temp_file_path)
    
    elif selected_option == "Extract sequence by random rate":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        random_rate = st.number_input("Random Rate", step=0.01)
        
        if file_path:
            if st.button("Run"):
                temp_file_path = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                out_file_path = random_rate_sequence(temp_file_path, random_rate)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name=f"{file_path.name.split('.')[0]}_random_{random_rate}.fasta",
                    mime="application/octet-stream"
                )
                os.remove(out_file_path)
                os.remove(temp_file_path)

    elif selected_option == "Extract sequence by random number":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        random_number = st.number_input("Number of Sequences", step=1)
        
        if file_path:
            if st.button("Run"):
                temp_file_path = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                out_file_path = random_number_sequence(temp_file_path, random_number)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name=f"{file_path.name.split('.')[0]}_random_{random_number}.fasta",
                    mime="application/octet-stream"
                )
                os.remove(out_file_path)
                os.remove(temp_file_path)

    elif selected_option == "Extract sequence by group":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        
        if fasta_file and info_file:
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
            
            if st.button("Run"):
                with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                    tmp_fasta.write(fasta_file.read())
                    tmp_fasta_path = tmp_fasta.name
                zip_filename = group_sequence(tmp_fasta_path, info_file, id_column, type_column)
                st.success("File processed successfully!")
                with open(zip_filename, 'rb') as f:
                    bytes_data = f.read()
                st.download_button(
                    label="Download Output",
                    data=bytes_data,
                    file_name="output.zip"
                )

    elif selected_option == "Extract sequence by gene site":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        
        if fasta_file and info_file:
            column_names = get_column_names(info_file, "xlsx")
            gene_column = st.selectbox("Gene Column Name", column_names)
            start_column = st.selectbox("Start Column Name", column_names)
            end_column = st.selectbox("End Column Name", column_names)
            
            if st.button("Run"):
                with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                    tmp_fasta.write(fasta_file.read())
                    tmp_fasta_path = tmp_fasta.name
                zip_filename = extract_gene(tmp_fasta_path, info_file, gene_column, start_column, end_column)
                st.success("File processed successfully!")
                with open(zip_filename, 'rb') as f:
                    bytes_data = f.read()
                st.download_button(
                    label="Download Output",
                    data=bytes_data,
                    file_name="output.zip")

    elif selected_option == "Split FASTA file into multiple files":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        sequences_per_file = st.number_input("Sequences per File", min_value=1, step=1)
        
        if file_path:
            if st.button("Run"):
                temp_file_path = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                zip_filename = split_fasta(temp_file_path, sequences_per_file)
                st.success("File processed successfully!")
                with open(zip_filename, 'rb') as f:
                    bytes_data = f.read()
                st.download_button(
                    label="Download Output",
                    data=bytes_data,
                    file_name="split_output.zip"
                )
                os.remove(temp_file_path)

if __name__ == "__main__":
    main()
