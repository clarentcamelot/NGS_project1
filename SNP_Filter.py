import argparse
from Bio import SeqIO

class FastaProcessor:
    def __init__(self, input_file, output_file, vcf_file):
        self.input_file = input_file
        self.output_file = output_file
        self.vcf_file = vcf_file
        self.variant_positions = set()
        self.filtered_sequences = []

    def process_vcf(self):
        with open(self.vcf_file, "r") as vcf_file:
            for line in vcf_file:
                if not line.startswith("#"): # Ignore comment lines
                    variant_position = line.split("\t")[0] # Extract variant position
                    self.variant_positions.add(variant_position)

    def filter_sequences(self):
        for record in SeqIO.parse(self.input_file, "fasta"):
            if record.id in self.variant_positions: # Check if sequence ID is in variant positions
                self.filtered_sequences.append(record)

    def write_filtered_fasta(self):
        with open(self.output_file, "w") as output_file:
            SeqIO.write(self.filtered_sequences, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Process FASTA file based on VCF information")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", dest="output_file", required=True, help="Output FASTA file")
    parser.add_argument("-v", "--vcf", dest="vcf_file", required=True, help="VCF file containing variant information")
    args = parser.parse_args()

    processor = FastaProcessor(args.input_file, args.output_file, args.vcf_file)
    processor.process_vcf()
    processor.filter_sequences()
    processor.write_filtered_fasta()

if __name__ == "__main__":
    main()
mkdir ../vc
