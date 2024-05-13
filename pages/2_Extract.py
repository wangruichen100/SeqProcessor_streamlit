import os
import io
import streamlit as st
import pandas as pd
import tempfile
import random
import shutil

# 导入要用于处理的函数
from utils import fasta_read, process_and_download2, get_column_names

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

# 随机抽取一定比例的序列
def random_rate_sequence(file_path, random_rate):
    sequences = fasta_read(file_path)
    # 计算需要选择的条目数量
    num_items_to_select = int(len(sequences) * random_rate)

    # 将字典视图转换为列表
    sequence_list = list(sequences.items())

    # 随机选择指定数量的条目
    sequences_random = random.sample(sequence_list, num_items_to_select)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, seq in sequences_random:
                f.write(f'>{name}\n{seq}\n')
    return out_path

# 随机抽取一定数量的序列
def random_number_sequence(file_path, num_sequences):
    sequences = fasta_read(file_path)

    # 获取所有的序列名
    sequence_names = list(sequences.keys())

    # 随机选择指定数量的序列名
    selected_names = random.sample(sequence_names, num_sequences)

    # 构建选定序列的字典
    selected_sequences = {name: sequences[name] for name in selected_names}

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as outfile:
        out_path = outfile.name
        with open(out_path, 'w') as f:
            for name, seq in selected_sequences.items():
                f.write(f'>{name}\n{seq}\n')
    return out_path

# 根据分组抽取序列
def group_sequence(file_path, info_path, id_column, type_column):
    seq_record = fasta_read(file_path) 
    df = pd.read_excel(info_path)
    seq_type = sorted(set(df[type_column].to_list()))
    output_folder = tempfile.mkdtemp()

    for seq_type_ in seq_type:
        df_type = df[df[type_column]==seq_type_]
        with open(os.path.join(output_folder, f"{seq_type_}.fas"), "w") as f:
            for name_ in df_type[id_column].to_list():          
                f.write(f">{name_}\n{seq_record[name_]}\n")

    zip_filename = shutil.make_archive(output_folder, 'zip', output_folder)
    shutil.rmtree(output_folder)
    return zip_filename


def main():
    st.title("Sequence Extraction")
# 创建顶部选项
    selected_option = st.radio("Select an option", ["Extract sequence by id",
                                                    "Extract sequence by random rate",
                                                    "Extract sequence by random number",
                                                    "Extract sequence by group"])
    
    if selected_option == "Extract sequence by id":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        info_path = st.file_uploader("Upload an Info file", type=["csv", "table", "txt", "xlsx"])
        
        process_and_download2(file_path, info_path, lambda x, y: extract_sequence(x, y), "require_seq.fasta")  # Define the output filename here
    
    elif selected_option == "Extract sequence by random rate":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        random_rate = st.number_input("Random Rate", step = 0.01)
        if file_path:
            if st.button("Run"):
                # Save uploaded file to a temporary location
                temp_file_path = tempfile.NamedTemporaryFile(delete=True, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                out_file_path = random_rate_sequence(temp_file_path, random_rate)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name= f"{file_path.name.split(".")[0]}_random_{random_rate}.fasta",  # Use the defined output filename here
                    mime="application/octet-stream"
                )
                os.remove(out_file_path)
                os.remove(temp_file_path)

    elif selected_option == "Extract sequence by random number":
        file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas", "fa"])
        random_number = st.number_input("Random Rate", step = 1)
        if file_path:
            if st.button("Run"):
                # Save uploaded file to a temporary location
                temp_file_path = tempfile.NamedTemporaryFile(delete=True, suffix=".fasta").name
                with open(temp_file_path, "wb") as f:
                    f.write(file_path.getvalue())
                
                out_file_path = random_number_sequence(temp_file_path, random_number)
                st.success("File processed successfully!")
                with open(out_file_path, "r") as f:
                    processed_content = f.read()
                st.download_button(
                    label="Download Output",
                    data=processed_content,
                    file_name= f"{file_path.name.split(".")[0]}_random_{random_number}.fasta",  # Use the defined output filename here
                    mime="application/octet-stream"
                )
                os.remove(out_file_path)
                os.remove(temp_file_path)

    elif selected_option == "Extract sequence by group":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fasta","fas"])
        info_file = st.file_uploader("Upload Info file (Excel)", type=["xlsx"])
        if fasta_file is not None and info_file is not None:
            # 获取上传的文件的列名
            column_names = get_column_names(info_file, "xlsx")
            id_column = st.selectbox("ID Column Name", column_names)
            type_column = st.selectbox("Type Column Name", column_names)
        if st.button("Run"):
            # 将上传的文件对象转换为文件路径
            with tempfile.NamedTemporaryFile(delete=False) as tmp_fasta:
                tmp_fasta.write(fasta_file.read())
                tmp_fasta_path = tmp_fasta.name
            zip_filename = group_sequence(tmp_fasta_path, info_file, id_column, type_column)
            with open(zip_filename, 'rb') as f:
                bytes_data = f.read()
                st.download_button(label="Download Output", data=bytes_data, file_name="output.zip")

if __name__ == "__main__":
    main()
