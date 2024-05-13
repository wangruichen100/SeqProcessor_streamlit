import os
import io
import streamlit as st
import pandas as pd
from Bio import SeqIO
from collections import Counter
from enum import Enum
import tempfile

# 导入要用于处理的函数
from utils import fasta_read, process_and_download2

# 提取序列
def extract_sequence(file_path, info_path):

    sequences = fasta_read(file_path)
    df_info = pd.read_csv(info_path)
    require_id = df_info.iloc[:, 0].to_list()  # 加上括号

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


def main():
    st.title("Sequence Extraction")
# 创建顶部选项
    selected_option = st.radio("Select an option", ["Extract sequence by id",
                                                    ])
    
    if selected_option == "Extract sequence by id":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        info_path = st.file_uploader("Upload an Info file", type=["csv", "table", "txt", "xlsx"])
        
        process_and_download2(file_path, info_path, lambda x, y: extract_sequence(x, y), "require_seq.fasta")  # Define the output filename here

if __name__ == "__main__":
    main()
