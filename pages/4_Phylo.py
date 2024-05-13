import subprocess
import streamlit as st
import os
import zipfile
from collections import deque
import tempfile
from utils import delete_folder_contents

# Check if IQ-TREE is installed
def is_iqtree_installed():
    try:
        # Execute conda list command and capture the output
        output = subprocess.check_output(["conda", "list"], universal_newlines=True)
        # Check if "IQ-TREE" is in the output
        return "iqtree" in output
    except subprocess.CalledProcessError:
        # If the command fails, Conda may not be installed or configured correctly
        return False

# Install IQ-TREE using Conda
def install_iqtree_with_conda():
    try:
        # Try installing IQ-TREE using the correct Conda installation command
        subprocess.run(["conda", "install", "bioconda::iqtree", "-y"], check=True)
        return True
    except subprocess.CalledProcessError:
        return False

# Command to perform phylogenetic tree reconstruction with IQ-TREE
def tree_iqtree(file_path: str, model: str, threads: str, output_placeholder):
    output_placeholder.empty()  # Clear the previous output

    output_text = ""
    output_text += "Check if IQ-TREE is installed\n"
    if not is_iqtree_installed():
        output_text += "IQ-TREE is not installed. Attempting to install with Conda...\n"
        if not install_iqtree_with_conda():
            output_text += "Failed to install IQ-TREE with Conda.\n"
            return
        output_text += "IQ-TREE installed successfully with Conda.\n"

    # Execute IQ-TREE command to build the tree
    # "iqtree -s input.fas -m MFP -bb 1000 -T AUTO"
    output_text += "Tree construction started...\n"
    try:
        process = subprocess.Popen(["iqtree", "-s", file_path, "-m", model, "-bb", "1000", "-T", threads],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        last_10_lines = deque(maxlen=10)  # Maintain the last 10 lines of output

        for line in process.stdout:
            last_10_lines.append(line)
            output_text += line
            output_placeholder.text("".join(last_10_lines))
        for line in process.stderr:
            last_10_lines.append(line)
            output_text += line
            output_placeholder.text("".join(last_10_lines))
        output_text += "Tree construction completed.\n"
    except subprocess.CalledProcessError as e:
        output_text += f"Error building tree with IQ-TREE: {e}\n"

    return output_text

def main():
    st.title("Phylogenetic Tree Reconstruction with IQ-TREE")

    file_path = st.file_uploader("Upload a FASTA file", type=["fasta", "fas"])
    if file_path is not None:
        model = st.selectbox("Select model", ["MFP", "GTR", "HKY", "JC", "K80"])
        threads = st.selectbox("Select number of threads", ["AUTO", "1", "2", "4", "8"])
        
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save uploaded file to a temporary location
            temp_file_path = os.path.join(temp_dir, file_path.name)
            with open(temp_file_path, "wb") as f:
                f.write(file_path.getvalue())

            if st.button("Build Tree"):
                output_placeholder = st.empty()  # Placeholder for real-time output
                output_text = tree_iqtree(temp_file_path, model, threads, output_placeholder)
                output_placeholder.text(output_text)

                # Compress output files into a zip archive
                output_files = [os.path.join(temp_dir, filename) for filename in os.listdir(temp_dir)]
                zip_file_path = os.path.join(temp_dir, f"{file_path.name}_iqtree.zip")
                with zipfile.ZipFile(zip_file_path, "w") as zipf:
                    for file in output_files:
                        zipf.write(file, os.path.basename(file))

                # Download button
                if os.path.exists(zip_file_path):
                    st.download_button(label="Download Output", data=open(zip_file_path, "rb"), file_name=zip_file_path.replace(temp_dir + "/", ""))

    # Remove the temporary file
    # delete_folder_contents(temp_dir)

if __name__ == "__main__":
    main()
