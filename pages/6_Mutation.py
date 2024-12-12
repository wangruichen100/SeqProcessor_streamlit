import streamlit as st
import pandas as pd
import os
import tempfile
import shutil
import zipfile

from Bio import SeqIO
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

from utils import fasta_read, fasta_read2, tab21_black, get_column_names

def pair_mutation(file_path, out_format):
    """
    Perform pairwise mutation analysis on sequences from a FASTA file and output results in specified format.

    Args:
        file_path (str): Path to the input FASTA file.
        out_format (str): Desired output format, either "csv", "table", or "excel".

    Returns:
        str: Path to the output ZIP file containing the results.
    """
    names, sequences = fasta_read2(file_path)
    ref_sequence = sequences[0]

    all_differences = []
    total_sequences = len(names) - 1  # Total sequences excluding the reference sequence

    # Create a progress bar widget
    progress_bar = st.progress(0)

    for idx, (name, sequence) in enumerate(zip(names[1:], sequences[1:]), 1):
        # Calculate progress percentage
        progress_percent = idx / total_sequences

        # Update the progress bar
        progress_bar.progress(int(progress_percent * 100))

        differences = find_differences(ref_sequence, sequence, name)
        all_differences.extend(differences)
    
    df = pd.DataFrame(all_differences, columns=["Site", "Ref seq", "Aligned seq", "Aligned seq name"])
    wide_df = df.pivot(index='Site', columns='Aligned seq name', values='Aligned seq')
    
    ref_sites = sorted(set(df["Site"]))
    ref_letters = [ref_sequence[i - 1] for i in ref_sites]
    wide_df.insert(0, 'Ref seq', ref_letters)

    for index, row in wide_df.iterrows():
        if pd.isna(row['Ref seq']):
            continue
        for col in wide_df.columns[1:]:
            if pd.isna(row[col]):
                wide_df.at[index, col] = row['Ref seq']
    
    tmp_dir = tempfile.mkdtemp()
    if out_format == "csv":
        df.to_csv(os.path.join(tmp_dir, "pair_mutation_l.csv"), index=False)
        wide_df.to_csv(os.path.join(tmp_dir, "pair_mutation_w.csv"))
    elif out_format == "table":
        df.to_csv(os.path.join(tmp_dir, "pair_mutation_l.tab"), sep="\t", index=False)
        wide_df.to_csv(os.path.join(tmp_dir, "pair_mutation_w.tab"), sep="\t")
    elif out_format == "excel":
        df.to_excel(os.path.join(tmp_dir, "pair_mutation_l.xlsx"), index=False)
        wide_df.to_excel(os.path.join(tmp_dir, "pair_mutation_w.xlsx"))

    zip_filename = shutil.make_archive(tmp_dir, 'zip', tmp_dir)
    shutil.rmtree(tmp_dir)

    return zip_filename

def find_differences(ref_seq, seq, name):
    """
    Find differences between reference sequence and given sequence.

    Args:
        ref_seq (str): Reference sequence.
        seq (str): Sequence to compare.
        name (str): Name of the sequence.

    Returns:
        list: List of differences as tuples (site, ref_base, seq_base, name).
    """
    differences = []
    for i, (ref_base, seq_base) in enumerate(zip(ref_seq, seq)):
        if ref_base != seq_base:
            differences.append([i+1, ref_base, seq_base, name])
    return differences

def group_aamutation(file_path, info_path, id_column, type_column, picture_type):
    """
    Perform group amino acid mutation analysis.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to the input info file.
        id_column (str): Column name for IDs.
        type_column (str): Column name for types.
        picture_type (str): Type of picture to generate ("bar" or "line").

    Returns:
        str: Path to the output ZIP file containing the results.
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

    pdf = os.path.join(output_folder, "group_aamutation.pdf")
    site_all, aa_all, stat_df_all, freq_df_all, group_all = [], [], [], [], []

    for name in seq_types:
        site, aa, stat_df, freq_df, group = aatj(os.path.join(output_folder, f"{name}.fas"), name)
        site_all += site
        aa_all += aa
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group

    site_all = [x + 1 for x in site_all]

    df = pd.DataFrame({
        "site": site_all,
        "aa": aa_all,
        "stat": stat_df_all,
        "freq": freq_df_all,
        "group": group_all
    }).fillna(0)
    df.to_csv(os.path.join(output_folder, "group_aamutation.csv"), index=False)

    with PdfPages(pdf) as pdf_pages:
        progress_bar = st.progress(0)
        for i in range(df["site"].min(), df["site"].max() + 1):
            df_site = df[df["site"] == i]
            pivot_df = df_site.pivot(index='group', columns='aa', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap=tab21_black(), yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap=tab21_black(), yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1), loc='upper right', fontsize=5)
            plt.title('Site' + str(i))
            plt.xlabel(None)
            pdf_pages.savefig()
            plt.clf()
            plt.close('all')

            progress = (i - df["site"].min() + 1) / (df["site"].max() - df["site"].min() + 1)
            progress_bar.progress(int(progress * 100))

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename

def aatj(input_fasta, group):
    """
    Perform amino acid type join analysis.

    Args:
        input_fasta (str): Path to the input FASTA file.
        group (str): Group name.

    Returns:
        tuple: Lists containing site, aa, stat_df, freq_df, and group.
    """
    site, aa, stat_df, freq_df = [], [], [], []

    seq_data = [list(record.seq) for record in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    aa_list = ['G', 'A', 'V', 'L', 'I', 'F', 'W', 'Y', 'D', 'N', 'E', 'K', 'Q', 'M', 'S', 'T', 'C', 'P', 'H', 'R', "-"]
    for column in df.columns:
        aa += aa_list
        site += [column] * len(aa_list)

        stat = Counter(df[column])
        aa_dict = {aa: stat.get(aa, 0) for aa in aa_list}
        stat_df.extend(aa_dict.values())

        freq = np.array(list(aa_dict.values())) / sum(aa_dict.values())
        freq_df.extend(freq)

    group_list = [group] * len(site)

    return site, aa, stat_df, freq_df, group_list

def group_ntmutation(file_path, info_path, id_column, type_column, picture_type):
    """
    Perform group nucleotide mutation analysis.

    Args:
        file_path (str): Path to the input FASTA file.
        info_path (str): Path to the input info file.
        id_column (str): Column name for IDs.
        type_column (str): Column name for types.
        picture_type (str): Type of picture to generate ("bar" or "line").

    Returns:
        str: Path to the output ZIP file containing the results.
    """
    seq_record = fasta_read(file_path)
    df = pd.read_excel(info_path)
    seq_types = sorted(set(df[type_column].to_list()))
    output_folder = tempfile.mkdtemp()

    for seq_type in seq_types:
        df_type = df[df[type_column] == seq_type]
        with open(os.path.join(output_folder, f"{seq_type}.fas"), "w") as f:
            for name in df_type[id_column].to_list():
                f.write(f">{name}\n{seq_record[name].upper()}\n")

    pdf = os.path.join(output_folder, "group_ntmutation.pdf")
    site_all, nt_all, stat_df_all, freq_df_all, group_all = [], [], [], [], []

    for name in seq_types:
        site, nt, stat_df, freq_df, group = nttj(os.path.join(output_folder, f"{name}.fas"), name)
        site_all += site
        nt_all += nt
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group

    site_all = [x + 1 for x in site_all]

    df = pd.DataFrame({
        "site": site_all,
        "nt": nt_all,
        "stat": stat_df_all,
        "freq": freq_df_all,
        "group": group_all
    }).fillna(0)
    df.to_csv(os.path.join(output_folder, "group_ntmutation.csv"), index=False)

    with PdfPages(pdf) as pdf_pages:
        progress_bar = st.progress(0)
        for i in range(df["site"].min(), df["site"].max() + 1):
            df_site = df[df["site"] == i]
            pivot_df = df_site.pivot(index='group', columns='nt', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap="Set1", yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap="Set1", yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1), loc='upper right', fontsize=5)
            plt.title('Site' + str(i))
            plt.xlabel(None)
            pdf_pages.savefig()
            plt.clf()
            plt.close('all')

            progress = (i - df["site"].min() + 1) / (df["site"].max() - df["site"].min() + 1)
            progress_bar.progress(int(progress * 100))

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename

def nttj(input_fasta, group):
    """
    Perform nucleotide type join analysis.

    Args:
        input_fasta (str): Path to the input FASTA file.
        group (str): Group name.

    Returns:
        tuple: Lists containing site, nt, stat_df, freq_df, and group.
    """
    site, nt, stat_df, freq_df = [], [], [], []

    seq_data = [list(record.seq) for record in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    nt_list = ['A', 'T', 'C', 'G', '-']
    for column in df.columns:
        nt += nt_list
        site += [column] * len(nt_list)

        stat = Counter(df[column])
        nt_dict = {nt: stat.get(nt, 0) for nt in nt_list}
        stat_df.extend(nt_dict.values())

        freq = np.array(list(nt_dict.values())) / sum(nt_dict.values())
        freq_df.extend(freq)

    group_list = [group] * len(site)

    return site, nt, stat_df, freq_df, group_list

def main():
    """
    Main function to run the Streamlit app for mutation analysis.
    """
    st.title("Mutation Analysis App")

    st.header("Analysis Options")
    analysis_type = st.radio("Select Analysis Type", ["Pair Mutation", "Group Amino Acid Mutation", "Group Nucleotide Mutation"])

    st.header("Upload Files")
    
    if analysis_type == "Pair Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas"])
        if fasta_file is not None:
            out_format = st.selectbox("Output Format", ["csv", "table", "excel"])

    elif analysis_type == "Group Amino Acid Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        if fasta_file is not None and info_file is not None:
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
            picture_type = st.selectbox("Picture Type", ["bar", "line"])

    elif analysis_type == "Group Nucleotide Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fas"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        if fasta_file is not None and info_file is not None:
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
            picture_type = st.selectbox("Picture Type", ["bar", "line"])

    if st.button("Run"):
        if analysis_type == "Pair Mutation":
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = pair_mutation(tmp_fasta_path, out_format)
            
        elif analysis_type == "Group Amino Acid Mutation":
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = group_aamutation(tmp_fasta_path, info_file, id_column, type_column, picture_type)
        
        elif analysis_type == "Group Nucleotide Mutation":
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = group_ntmutation(tmp_fasta_path, info_file, id_column, type_column, picture_type)

        st.success("Analysis completed successfully. Please download the output zip file.")
        with open(zip_filename, 'rb') as f:
            bytes_data = f.read()
            st.download_button(label="Download Output", data=bytes_data, file_name="output.zip")

# Run the Streamlit app
if __name__ == "__main__":
    main()
