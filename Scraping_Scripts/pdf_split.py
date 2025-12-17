import os
from PyPDF2 import PdfReader, PdfWriter

def split_pdf(input_pdf, output_folder, pages_per_split=4):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Open the PDF
    reader = PdfReader(input_pdf)
    total_pages = len(reader.pages)
    
    # Split into chunks
    for i in range(0, total_pages, pages_per_split):
        writer = PdfWriter()
        for j in range(pages_per_split):
            if i + j < total_pages:
                writer.add_page(reader.pages[i + j])
        
        # Output file name
        output_filename = os.path.join(output_folder, f"split_{i//pages_per_split + 1}.pdf")
        with open(output_filename, "wb") as output_pdf:
            writer.write(output_pdf)
        print(f"Saved: {output_filename}")

# Example usage
split_pdf("input.pdf", "output_splits")
