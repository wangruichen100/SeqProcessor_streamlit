import os
import io
import streamlit as st
import pandas as pd
from collections import Counter
from enum import Enum
import tempfile
from collections import defaultdict

from utils import fasta_read2

def analyze_codon_sequence(sequence):
    sequence = sequence.upper().replace('T', 'U')

    sequence_length = len(sequence)
    
    if sequence_length % 3 != 0:
        raise ValueError("Sequence length is not a multiple of 3")
    
    codons = [sequence[i:i+3] for i in range(0, sequence_length, 3)]
    codons = [codon for codon in codons if all(base in 'AUCG' for base in codon)]
    
    codon_count = len(codons)
    
    A_count = sequence.count('A')
    U_count = sequence.count('U')
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    
    total_bases = A_count + U_count + C_count + G_count

    A_ratio = A_count / total_bases
    U_ratio = U_count / total_bases
    C_ratio = C_count / total_bases
    G_ratio = G_count / total_bases

    AU_ratio = A_ratio + U_ratio
    GC_ratio = G_ratio + C_ratio

    A1_ratio = sum(1 for codon in codons if codon[0] == 'A') / codon_count
    U1_ratio = sum(1 for codon in codons if codon[0] == 'U') / codon_count
    C1_ratio = sum(1 for codon in codons if codon[0] == 'C') / codon_count
    G1_ratio = sum(1 for codon in codons if codon[0] == 'G') / codon_count

    A2_ratio = sum(1 for codon in codons if codon[1] == 'A') / codon_count
    U2_ratio = sum(1 for codon in codons if codon[1] == 'U') / codon_count
    C2_ratio = sum(1 for codon in codons if codon[1] == 'C') / codon_count
    G2_ratio = sum(1 for codon in codons if codon[1] == 'G') / codon_count

    A3_ratio = sum(1 for codon in codons if codon[2] == 'A') / codon_count
    U3_ratio = sum(1 for codon in codons if codon[2] == 'U') / codon_count
    C3_ratio = sum(1 for codon in codons if codon[2] == 'C') / codon_count
    G3_ratio = sum(1 for codon in codons if codon[2] == 'G') / codon_count

    AU1_ratio = A1_ratio + U1_ratio
    GC1_ratio = G1_ratio + C1_ratio
    AU2_ratio = A2_ratio + U2_ratio
    GC2_ratio = G2_ratio + C2_ratio
    AU3_ratio = A3_ratio + U3_ratio
    GC3_ratio = G3_ratio + C3_ratio
    AU12_ratio = (AU1_ratio + AU2_ratio) / 2
    GC12_ratio = (GC1_ratio + GC2_ratio) / 2
    A3_AU3 = A3_ratio / AU3_ratio if AU3_ratio != 0 else 0
    G3_GC3 = G3_ratio / GC3_ratio if GC3_ratio != 0 else 0
    
    return {
        "Sequence length": sequence_length,
        'Effective nucleotide quantity': total_bases,
        'Effective codon quantity': codon_count,
        'A_ratio': A_ratio,
        'U_ratio': U_ratio,
        'C_ratio': C_ratio,
        'G_ratio': G_ratio,
        'AU_ratio': AU_ratio,
        'GC_ratio': GC_ratio,
        'A1_ratio': A1_ratio,
        'U1_ratio': U1_ratio,
        'C1_ratio': C1_ratio,
        'G1_ratio': G1_ratio,
        'A2_ratio': A2_ratio,
        'U2_ratio': U2_ratio,
        'C2_ratio': C2_ratio,
        'G2_ratio': G2_ratio,
        'A3_ratio': A3_ratio,
        'U3_ratio': U3_ratio,
        'C3_ratio': C3_ratio,
        'G3_ratio': G3_ratio,
        'AU1_ratio': AU1_ratio,
        'GC1_ratio': GC1_ratio,
        'AU2_ratio': AU2_ratio,
        'GC2_ratio': GC2_ratio,
        'AU3_ratio': AU3_ratio,
        'GC3_ratio': GC3_ratio,
        'AU12_ratio': AU12_ratio,
        'GC12_ratio': GC12_ratio,
        'A3_AU3': A3_AU3,
        'G3_GC3': G3_GC3
    }

def base_compsition(file_path):
    names, sequences = fasta_read2(file_path)
    df_values = []
    for name_, sequence_ in zip(names, sequences):
        df_values.append([name_] + list(analyze_codon_sequence(sequence_).values()))

    df = pd.DataFrame(df_values, columns=[
        "Name",
        "Sequence length",
        'Effective nucleotide quantity',
        'Effective codon quantity',
        'A_ratio',
        'U_ratio',
        'C_ratio',
        'G_ratio',
        'AU_ratio',
        'GC_ratio',
        'A1_ratio',
        'U1_ratio',
        'C1_ratio',
        'G1_ratio',
        'A2_ratio',
        'U2_ratio',
        'C2_ratio',
        'G2_ratio',
        'A3_ratio',
        'U3_ratio',
        'C3_ratio',
        'G3_ratio',
        'AU1_ratio',
        'GC1_ratio',
        'AU2_ratio',
        'GC2_ratio',
        'AU3_ratio',
        'GC3_ratio',
        'AU12_ratio',
        'GC12_ratio',
        'A3_AU3',
        'G3_GC3'
    ])
    
    return df



# 定义同义密码子表
codon_table = {
    'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Y': ['UAU', 'UAC'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['UGU', 'UGC'],
    'W': ['UGG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG']
}

# 计算RSCU值
def calculate_rscu(sequence):
    sequence = sequence.upper().replace('T', 'U')
    # 初始化同义密码子的计数字典
    codon_count = defaultdict(int)
    aa_codon_count = defaultdict(lambda: defaultdict(int))
    
    # 统计每个密码子出现的频率
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            for aa, codons in codon_table.items():
                if codon in codons:
                    codon_count[codon] += 1
                    aa_codon_count[aa][codon] += 1
                    break
    
    # 计算RSCU值
    rscu_values = defaultdict(float)
    for aa, codons in aa_codon_count.items():
        total_count = sum(codons.values())
        num_codons = len(codon_table[aa])
        for codon, count in codons.items():
            if total_count > 0:
                rscu_values[codon] = (count * num_codons) / total_count
            else:
                rscu_values[codon] = 0.0
    
    return rscu_values

def rscu(file_path):
    names, sequences = fasta_read2(file_path)
    all_rscu_values = []

    for name, sequence in zip(names, sequences):
        rscu_values = calculate_rscu(sequence)
        rscu_values['Name'] = name
        all_rscu_values.append(rscu_values)
    
    order = ["Name","UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG","UAU","UAC","UGU","UGC","UGG",
         "CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG","CAU","CAC","CAA","CAG","CGU","CGC","CGA",
         "CGG","AUU","AUC","AUA","AUG","ACU","ACC","ACA","ACG","AAU","AAC","AAA","AAG","AGU","AGC",
         "AGA","AGG","GUU","GUC","GUA","GUG","GCU","GCC","GCA","GCG","GAU","GAC","GAA","GAG","GGU",
         "GGC","GGA","GGG"]

    df = pd.DataFrame(all_rscu_values)
    df = df.reindex(columns=order)
    return df

def main():
    st.title("Codon Usage")
    selected_option = st.radio("Select an option", ["Base composition", "Relative Synonymous Codon Usage"])
    
    file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
    
    if file_path is not None:
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_file.write(file_path.read())
            temp_file_path = temp_file.name

        if st.button("Run"):
            if selected_option == "Base composition":
                df = base_compsition(temp_file_path)
                st.session_state['df'] = df
                st.session_state['output_filename_prefix'] = "base_composition"
            elif selected_option == "Relative Synonymous Codon Usage":
                df = rscu(temp_file_path)
                st.session_state['df'] = df
                st.session_state['output_filename_prefix'] = "rscu"

            st.dataframe(df)

    if 'df' in st.session_state and 'output_filename_prefix' in st.session_state:
        download_format = st.selectbox("Select download format", ["CSV", "Excel", "Table"])
        output_filename_prefix = st.session_state['output_filename_prefix']

        if download_format == "CSV":
            csv = st.session_state['df'].to_csv(index=False)
            st.download_button(
                label="Download as CSV",
                data=csv,
                file_name=f"{output_filename_prefix}.csv",
                mime="text/csv"
            )

        elif download_format == "Table":
            tab_delimited = st.session_state['df'].to_csv(index=False, sep='\t')
            st.download_button(
                label="Download as Table",
                data=tab_delimited,
                file_name=f"{output_filename_prefix}.tsv",
                mime="text/tab-separated-values"
            )

        elif download_format == "Excel":
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                st.session_state['df'].to_excel(writer, index=False, sheet_name='CodonUsage')
            excel_buffer.seek(0)
            st.download_button(
                label="Download as Excel",
                data=excel_buffer,
                file_name=f"{output_filename_prefix}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

if __name__ == "__main__":
    main()
