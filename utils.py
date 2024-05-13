import os
import tempfile
import streamlit as st
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

def fasta_read(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_name = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    sequences[sequence_name] = sequence
                sequence_name = line[1:]
                sequence = ''
            else:
                sequence += line
        # Adding the last sequence
        if sequence_name is not None:
            sequences[sequence_name] = sequence

    return sequences

def fasta_read2(file_path):
    names = []
    sequences = []
    current_name = None
    current_sequence = ""
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    names.append(current_name)
                    sequences.append(current_sequence)
                current_name = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_name:
            names.append(current_name)
            sequences.append(current_sequence)
    return names, sequences

def delete_folder_contents(folder_path):
    try:
        # 遍历文件夹中的所有内容
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            # 如果是文件，则直接删除
            if os.path.isfile(item_path):
                os.remove(item_path)
            # 如果是文件夹，则递归调用删除子文件夹中的内容
            elif os.path.isdir(item_path):
                delete_folder_contents(item_path)
        print(f"Contents of folder {folder_path} deleted successfully.")
    except Exception as e:
        print(f"Error deleting contents of folder {folder_path}: {e}")

def translate_dna_to_protein(dna_sequence, codon_table):
    # Create a Seq object from the DNA sequence
    dna_seq = Seq(dna_sequence)
    
    # Translate the DNA sequence to protein sequence using the specified codon table
    protein_seq = dna_seq.translate(table=codon_table)
    
    return protein_seq

def translate_sequence(nucleotide_sequence):
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    
    nucleotide_sequence = nucleotide_sequence.upper()
    protein_sequence = ""
    for i in range(0, len(nucleotide_sequence), 3):
        codon = nucleotide_sequence[i:i+3]
        if codon in codon_table:
            protein_sequence += codon_table[codon]
        else:
            protein_sequence += "?"
    return protein_sequence


def tab21_black():
    # Get tab20 colormap
    tab20 = cm.get_cmap('tab20')

    tab20_black = tab20(np.linspace(0, 1, 20))
    tab20_black = np.vstack(([0, 0, 0, 1],tab20_black))  # [0, 0, 0] represents black color
    # Create colormap with black color added
    colormap = ListedColormap(tab20_black)
    return colormap

def dark2_black():
    # Get Dark2 colormap
    dark2 = cm.get_cmap('Dark2')

    dark2_black = dark2(np.linspace(0, 1, 8))
    dark2_black = np.vstack(([0, 0, 0, 1],dark2_black))  # [0, 0, 0] represents black color
    # Create colormap with black color added
    colormap = ListedColormap(dark2_black)
    return colormap

def get_column_names(file_path, file_format):
    if file_format == "xlsx":
        df = pd.read_excel(file_path, engine='openpyxl')
    elif file_format == "csv":
        df = pd.read_csv(file_path)
    elif file_format == "table" or file_format == "txt":
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError("Unsupported file format")

    return df.columns.tolist()

def process_and_download(file_path, process_function, output_filename):
    if file_path:
        if st.button("Run"):
            # Save uploaded file to a temporary location
            temp_file_path = tempfile.NamedTemporaryFile(delete=True, suffix=".fasta").name
            with open(temp_file_path, "wb") as f:
                f.write(file_path.getvalue())
            
            out_file_path = process_function(temp_file_path)
            st.success("File processed successfully!")
            with open(out_file_path, "r") as f:
                processed_content = f.read()
            st.download_button(
                label="Download Output",
                data=processed_content,
                file_name=output_filename,  # Use the defined output filename here
                mime="application/octet-stream"
            )
            os.remove(out_file_path)
            os.remove(temp_file_path)

def process_and_download2(file_path, info_path, process_function, output_filename):
    if file_path and info_path:
        if st.button("Run"):
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

            out_file_path = process_function(temp_fas_file_path, temp_info_file_path)
            st.success("File processed successfully!")
            with open(out_file_path, "r") as f:
                processed_content = f.read()
            st.download_button(
                label="Download Output",
                data=processed_content,
                file_name=output_filename,  # Use the defined output filename here
                mime="application/octet-stream"
            )

