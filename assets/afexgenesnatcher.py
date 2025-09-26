import os
import time
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# === SETTINGS ===
gff_file = "F:/AFEXBioLab/Database/AnthuriumAmnicola.final.gff/AnthuriumAmnicola.final.gff"
genome_file = "F:/AFEXBioLab/Database/AnthuriumAmnicola.genome.fa/AnthuriumAmnicola.genome.fa"
output_folder = "Extracted_Genes_GFF"
report_file = "lab_report.txt"

# === START ===
start_time = time.time()
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Load genome
genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Extract features from GFF
in_handle = open(gff_file)
out_records = []
saved_files = []
skipped = 0

for rec in GFF.parse(in_handle, base_dict=genome):
    for feature in rec.features:
        if feature.type in ["gene", "mRNA", "CDS"]:
            try:
                qualifiers = feature.qualifiers
                gene_id = qualifiers.get("ID", ["unnamed_gene"])[0].replace(":", "_")
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                seq = rec.seq[start:end]
                if strand == -1:
                    seq = seq.reverse_complement()

                gene_record = SeqRecord(
                    seq=seq,
                    id=gene_id,
                    description=f"{gene_id} | {start}-{end} | Strand: {strand}"
                )

                filename = f"{gene_id}.fasta"
                filepath = os.path.join(output_folder, filename)
                SeqIO.write(gene_record, filepath, "fasta")
                saved_files.append(filename)

            except Exception as e:
                skipped += 1
                print(f"❌ Skipped: {e}")

in_handle.close()
elapsed = round(time.time() - start_time, 2)

# === REPORT ===
with open(os.path.join(output_folder, report_file), "w", encoding="utf-8") as report:
    report.write("=== LAB REPORT 🧬 (from GFF) ===\n")
    report.write(f"✅ Total Genes Extracted: {len(saved_files)}\n")
    report.write(f"❌ Genes Skipped: {skipped}\n")
    report.write(f"📂 Output Folder: {output_folder}\n")
    report.write(f"📜 File Format: FASTA\n")
    if saved_files:
        report.write(f"🧬 Sample File: {saved_files[0]} ... {saved_files[-1]}\n")
    report.write(f"⏱️ Time Taken: {elapsed}s\n\n")
    report.write("🧪 Extracted Gene Files:\n")
    report.write("\n".join(saved_files))

print(f"\n🔥 DONE! Extracted {len(saved_files)} gene sequences from GFF + Genome!")
print(f"📄 LAB REPORT saved as '{report_file}'")
