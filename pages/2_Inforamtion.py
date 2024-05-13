import pandas as pd
import streamlit as st
import io

from Bio import SeqIO

def process_genebank_record(record):
    feature = record.features[0] if record.features else None
    qualifiers = {
        "id": record.name,
        "length": len(record.seq),
        "organism": feature.qualifiers.get("organism", ["no"])[0],
        "strain": feature.qualifiers.get("strain", ["no"])[0],
        "isolate": feature.qualifiers.get("isolate", ["no"])[0],
        "country": feature.qualifiers.get("country", ["no"])[0],
        "host": feature.qualifiers.get("host", ["no"])[0],
        "date": feature.qualifiers.get("collection_date", ["no"])[0],
        "note": feature.qualifiers.get("note", ["no"])[0],
        "taxonomy": record.annotations.get("taxonomy", "no"),
    }
    return qualifiers

def genbank_profile(file_path):
    with io.StringIO(file_path.getvalue().decode()) as handle:
        records = [process_genebank_record(record) for record in SeqIO.parse(handle, "genbank")]
    
    df = pd.DataFrame(records)
    
    return df

def main():
    st.title("Extracting information from Genbank file")
    file_path = st.file_uploader("Upload a Genbank file", type=["gb", "gbk"])
    
    out_format = st.selectbox("Output format", ["csv", "tab", "excel"])
    if out_format == "csv":
        out_extension = "csv"
    elif out_format == "tab":
        out_extension = "tab"
    elif out_format == "excel":
        out_extension = "xlsx"

    if file_path is not None:
        if st.button("Process"):
            df = genbank_profile(file_path)
            st.write(df)
            st.success("File processed successfully!")
        
            df = genbank_profile(file_path)  # Re-generate DataFrame
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
