import streamlit as st
import pandas as pd
import os
import tempfile  # 导入tempfile模块
import shutil
import zipfile

from Bio import SeqIO
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

from utils import fasta_read, fasta_read2, tab21_black

def get_column_names(file_path, file_format):
    if file_format == "xlsx":
        df = pd.read_excel(file_path)
    elif file_format == "csv":
        df = pd.read_csv(file_path)
    elif file_format == "table":
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError("Unsupported file format")

    return df.columns.tolist()

def pair_mutation(file_path, out_format):
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

    # Close the progress bar
    progress_bar.empty()

    return zip_filename

def find_differences(ref_seq, seq, name):
    differences = []
    for i, (letter1, letter2) in enumerate(zip(ref_seq, seq)):
        if letter1 != letter2:
            differences.append([i+1, letter1, letter2, name])
    return differences

def group_aamutation(file_path, info_path, id_column, type_column, picture_type):
    seq_record = fasta_read(file_path) 
    df = pd.read_excel(info_path)
    seq_type = sorted(set(df[type_column].to_list()))
    output_folder = tempfile.mkdtemp()

    for seq_type_ in seq_type:
        df_type = df[df[type_column]==seq_type_]
        with open(os.path.join(output_folder, f"{seq_type_}.fas"), "w") as f:
            for name_ in df_type[id_column].to_list():          
                f.write(f">{name_}\n{seq_record[name_]}\n")

    pdf = os.path.join(output_folder, "group_aamutation.pdf")

    site_all = []
    aa_all = []
    stat_df_all = []
    freq_df_all = []
    group_all = []

    for name_ in seq_type:
        site, aa, stat_df, freq_df, group_ = aatj(os.path.join(output_folder, f"{name_}.fas"), name_)

        site_all += site
        aa_all += aa
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group_

    df = pd.DataFrame()
    df["site"] = site_all
    df["aa"] = aa_all
    df["stat"] = stat_df_all
    df["freq"] = freq_df_all
    df["group"] = group_all
    df = df.fillna(0)
    df.to_csv(os.path.join(output_folder, "group_aamutation.csv"),index = False)

    with PdfPages(pdf) as pdf:
        progress_bar = st.progress(0)
        for i_ in range(df["site"].min(), df["site"].max() + 1):
            df_site = df[df["site"] == i_]

            pivot_df = df_site.pivot(index='group', columns='aa', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap=tab21_black(), yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap=tab21_black(), yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1),loc='upper right', fontsize=5)
            plt.title('Site' + str(i_+1))
            plt.xlabel(None)
            pdf.savefig()
            plt.clf()
            plt.close('all')

            # 更新进度条
            progress = (i_ - df["site"].min() + 1) / (df["site"].max() - df["site"].min() + 1)
            progress_bar.progress(int(progress * 100))

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename

def aatj(input_fasta, group):
    site = []
    aa = []
    stat_df = []
    freq_df = []

    seq_data = [list(i.seq) for i in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    for column in df.columns:
        aa_ = ['G', 'A', 'V', 'L', 'I', 'F', 'W', 'Y', 'D', 'N', 'E', 'K', 'Q', 'M', 'S', 'T', 'C', 'P', 'H', 'R', "-"]
        aa += aa_
        site += [column] * len(aa_)

        stat = Counter(df[column])
        aa_dict = {aa: stat.get(aa, 0) for aa in aa_}
        stat_df.extend(aa_dict.values())

        freq = np.array(list(aa_dict.values())) / sum(aa_dict.values())
        freq_df.extend(freq)

    group_ = [group] * len(site)

    return site, aa, stat_df, freq_df, group_

def group_ntmutation(file_path, info_path, id_column, type_column, picture_type):
    seq_record = fasta_read(file_path)
    df = pd.read_excel(info_path)
    seq_type = sorted(set(df[type_column].to_list()))
    output_folder = tempfile.mkdtemp()

    for seq_type_ in seq_type:
        df_type = df[df[type_column]==seq_type_]
        with open(os.path.join(output_folder, f"{seq_type_}.fas"), "w") as f:
            for name_ in df_type[id_column].to_list():          
                f.write(f">{name_}\n{seq_record[name_].upper()}\n")

    pdf = os.path.join(output_folder, "group_ntmutation.pdf")

    site_all = []
    nt_all = []
    stat_df_all = []
    freq_df_all = []
    group_all = []

    for name_ in seq_type:
        site, nt, stat_df, freq_df, group_ = nttj(os.path.join(output_folder, f"{name_}.fas"), name_)

        site_all += site
        nt_all += nt
        stat_df_all += stat_df
        freq_df_all += freq_df
        group_all += group_

    df = pd.DataFrame()
    df["site"] = site_all
    df["nt"] = nt_all
    df["stat"] = stat_df_all
    df["freq"] = freq_df_all
    df["group"] = group_all
    df = df.fillna(0)
    df.to_csv(os.path.join(output_folder, "group_ntmutation.csv"),index = False)

    with PdfPages(pdf) as pdf:
        progress_bar = st.progress(0)
        for i_ in range(df["site"].min(), df["site"].max() + 1):
            df_site = df[df["site"] == i_]

            pivot_df = df_site.pivot(index='group', columns='nt', values='freq')
            if picture_type == "bar":
                pivot_df.plot(kind='bar', stacked=True, colormap="Set1", yticks=[0.5, 1])
            elif picture_type == "line":
                pivot_df.plot(kind='line', colormap="Set1", yticks=[0.5, 1])
            plt.ylim((0, 1.1))
            plt.legend(prop={'size': 8}, bbox_to_anchor=(1.12, 1),loc='upper right', fontsize=5)
            plt.title('Site' + str(i_+1))
            plt.xlabel(None)
            pdf.savefig()
            plt.clf()
            plt.close('all')

            # 更新进度条
            progress = (i_ - df["site"].min() + 1) / (df["site"].max() - df["site"].min() + 1)
            progress_bar.progress(int(progress * 100))

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename

def nttj(input_fasta, group):
    site = []
    nt = []
    stat_df = []
    freq_df = []

    seq_data = [list(i.seq) for i in SeqIO.parse(input_fasta, "fasta")]
    df = pd.DataFrame(seq_data)

    for column in df.columns:
        nt_ = ['A',"T","C","G","-"]
        nt += nt_
        site += [column] * len(nt_)

        stat = Counter(df[column])
        nt_dict = {nt: stat.get(nt, 0) for nt in nt_}
        stat_df.extend(nt_dict.values())

        freq = np.array(list(nt_dict.values())) / sum(nt_dict.values())
        freq_df.extend(freq)

    group_ = [group] * len(site)

    return site, nt, stat_df, freq_df, group_

def main():
    st.title("Mutation Analysis App")

    st.header("Analysis Options")
    analysis_type = st.radio("Select Analysis Type", ["Pair Mutation", "Group Amino Acid Mutation", "Group Nucleotide Mutation"])

    st.header("Upload Files")
    
    if analysis_type == "Pair Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta","fas"])
        if fasta_file is not None:
            out_format = st.selectbox("Output Format", ["csv", "table", "excel"])

    elif analysis_type == "Group Amino Acid Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta","fas"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        if fasta_file is not None and info_file is not None:
            # 获取上传的文件的列名
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
            picture_type = st.selectbox("Picture Type", ["bar", "line"])

    elif analysis_type == "Group Nucleotide Mutation":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta","fas"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        if fasta_file is not None and info_file is not None:
        # 获取上传的文件的列名
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
            picture_type = st.selectbox("Picture Type", ["bar", "line"])

    if st.button("Run"):
        if analysis_type == "Pair Mutation":
            # 将上传的文件对象转换为文件路径
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = pair_mutation(tmp_fasta_path, out_format)
            
        elif analysis_type == "Group Amino Acid Mutation":
            # 将上传的文件对象转换为文件路径
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = group_aamutation(tmp_fasta_path, info_file, id_column, type_column, picture_type)
        elif analysis_type == "Group Nucleotide Mutation":
            # 将上传的文件对象转换为文件路径
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = group_ntmutation(tmp_fasta_path, info_file, id_column, type_column, picture_type)

        # 提示处理完成
        st.success("Analysis completed successfully. Please download the output zip file.")
        with open(zip_filename, 'rb') as f:
            bytes_data = f.read()
            st.download_button(label="Download Output", data=bytes_data, file_name="output.zip")

# Run the Streamlit app
if __name__ == "__main__":
    main()
