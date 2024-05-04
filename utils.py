import os


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

# 调用函数并传入要删除的文件夹路径
delete_folder_contents('/path/to/folder')
