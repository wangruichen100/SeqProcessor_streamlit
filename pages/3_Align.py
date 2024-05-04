import subprocess
import streamlit as st
import os
from collections import deque

# Command to perform multiple sequence alignment with MAFFT
# Check if MAFFT is installed
def is_mafft_installed():
    try:
        # Execute conda list command and capture the output
        output = subprocess.check_output(["conda", "list"], universal_newlines=True)
        # Check if "mafft" is in the output
        return "mafft" in output
    except subprocess.CalledProcessError:
        # If the command fails, Conda may not be installed or configured correctly
        return False

# Install MAFFT using Conda
def install_mafft_with_conda():
    try:
        # Try installing MAFFT using the correct Conda installation command
        subprocess.run(["conda", "install", "bioconda::mafft", "-y"], check=True)
        return True
    except subprocess.CalledProcessError:
        return False


def align_mafft(file_path, out_path, output_placeholder):
    output_placeholder.empty()  # Clear the previous output

    output_placeholder.text("Check if MAFFT is installed\n")
    if not is_mafft_installed():
        output_placeholder.text("MAFFT is not installed. Attempting to install with Conda...\n")
        if not install_mafft_with_conda():
            output_placeholder.text("Failed to install MAFFT with Conda.\n")
            return
        output_placeholder.text("MAFFT installed successfully with Conda.\n")

    # Execute the MAFFT command
    output_placeholder.text("Running MAFFT command...\n")
    try:
        # Execute the MAFFT command and capture the output
        process = subprocess.Popen(["mafft", file_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        last_10_lines = deque(maxlen=10)  # Maintain the last 10 lines of output
        output_text = ""
        for line in process.stdout:
            last_10_lines.append(line)
            output_text += line
            output_placeholder.text("".join(last_10_lines))

        # Write the output to the file in the output directory
        # with open(os.path.join(os.path.dirname(out_path), "mafft_output.txt"), "w") as f:
        #     f.write(output_text)
    except subprocess.CalledProcessError as e:
        output_placeholder.text(f"Error running MAFFT: {e}\n")
        return

    # Write the result to the output file
    # output_placeholder.text(f"Writing output to {out_path}\n")
    try:
        with open(out_path, "w") as f:
            subprocess.run(["mafft", file_path], check=True, stdout=f)
        output_placeholder.text("MAFFT run successful.\n")
    except subprocess.CalledProcessError as e:
        output_placeholder.text(f"Error writing output file: {e}\n")

def main():
    st.title("Multiple sequence alignment by MAFFT")
    file_path = st.file_uploader("Upload a Fasta file", type=["fasta", "fas"])
    if file_path is not None:
        if st.button("Process"):
            # Save uploaded file to a temporary location
            temp_file_path = "./temp/temp_file.fasta"
            with open(temp_file_path, "wb") as f:
                f.write(file_path.getvalue())

            out_path = f"./temp/{file_path.name.split('.')[0]}_mafft.fas"
            output_placeholder = st.empty()  # Placeholder for real-time output
            align_mafft(temp_file_path, out_path, output_placeholder)
            st.success("File processed successfully!")  

            # Remove the temporary file
            os.remove(temp_file_path)

            # Download button
            if os.path.exists(out_path):
                st.download_button(label="Download Result", data=open(out_path, "rb"), file_name=out_path.replace("./temp/",""))
                os.remove(out_path)

if __name__ == "__main__":
    main()
