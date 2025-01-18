import pandas as pd
import streamlit as st
import io
from Bio import SeqIO

def process_genebank_record(record):
    """
    Process a single GenBank record to extract relevant information.

    Args:
        record (SeqRecord): A Biopython SeqRecord object from a GenBank file.

    Returns:
        dict: A dictionary containing extracted qualifiers and annotations.
    """
    feature = record.features[0] if record.features else None
    qualifiers = {
        "id": record.name,
        "length": len(record.seq),
        "organism": feature.qualifiers.get("organism", ["no"])[0] if feature else "no",
        "strain": feature.qualifiers.get("strain", ["no"])[0] if feature else "no",
        "isolate": feature.qualifiers.get("isolate", ["no"])[0] if feature else "no",
        "country": feature.qualifiers.get("geo_loc_name", ["no"])[0] if feature else "no",
        "host": feature.qualifiers.get("host", ["no"])[0] if feature else "no",
        "date": feature.qualifiers.get("collection_date", ["no"])[0] if feature else "no",
        "note": feature.qualifiers.get("note", ["no"])[0] if feature else "no",
        "taxonomy": record.annotations.get("taxonomy", "no"),
    }
    return qualifiers

def genbank_profile(file_path):
    """
    Generate a DataFrame from a GenBank file.

    Args:
        file_path (UploadedFile): The uploaded GenBank file.

    Returns:
        DataFrame: A DataFrame containing extracted information from the GenBank file.
    """
    with io.StringIO(file_path.getvalue().decode()) as handle:
        records = [process_genebank_record(record) for record in SeqIO.parse(handle, "genbank")]
    
    df = pd.DataFrame(records)
    return df

def main():
    """
    Main function to run the Streamlit app for extracting information from GenBank files.
    """
    st.title("Extracting Information from GenBank File")
    file_path = st.file_uploader("Upload a GenBank file", type=["gb", "gbk"])
    
    out_format = st.selectbox("Output Format", ["csv", "tab", "excel"])
    out_extension = {
        "csv": "csv",
        "tab": "tab",
        "excel": "xlsx"
    }[out_format]

    if file_path is not None:
        if st.button("Run"):
            df = genbank_profile(file_path)
            st.write(df)
            st.success("File processed successfully!")
        
            output_filename = file_path.name.split('.')[0] + f'.{out_extension}'
            output = io.BytesIO()
            if out_format == "csv":
                df.to_csv(output, index=False)
            elif out_format == "tab":
                df.to_csv(output, sep="\t", index=False)
            elif out_format == "excel":
                df.to_excel(output, index=False)
            output.seek(0)
            st.download_button(
                label="Download Output",
                data=output,
                file_name=output_filename,
                mime="application/octet-stream"
            )

if __name__ == "__main__":
    main()
