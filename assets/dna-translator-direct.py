#!/usr/bin/env python3
"""
AFEX Genesis™ - Direct DNA Translator (No API)
Command-line DNA translation using BioPython
"""

import sys
import json
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import io

class DNATranslatorDirect:
    """Direct DNA translation without API"""
    
    def __init__(self):
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def validate_dna_sequence(self, sequence):
        """Validate DNA sequence contains only valid nucleotides"""
        valid_bases = set('ATGCN')
        sequence_upper = sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '').replace('\t', '')
        return all(base in valid_bases for base in sequence_upper), sequence_upper
    
    def parse_file_content(self, file_content, file_format):
        """Parse different file formats and extract DNA sequences"""
        sequences = []
        
        try:
            if file_format.lower() in ['fasta', 'fa']:
                # Parse FASTA format
                fasta_io = io.StringIO(file_content)
                for record in SeqIO.parse(fasta_io, "fasta"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq)
                    })
            
            elif file_format.lower() in ['genbank', 'gb']:
                # Parse GenBank format
                genbank_io = io.StringIO(file_content)
                for record in SeqIO.parse(genbank_io, "genbank"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq)
                    })
            
            elif file_format.lower() == 'txt':
                # Parse plain text - assume it's just DNA sequence
                clean_sequence = file_content.replace(' ', '').replace('\n', '').replace('\r', '').replace('\t', '')
                if clean_sequence:
                    sequences.append({
                        'id': 'sequence_1',
                        'description': 'Plain text sequence',
                        'sequence': clean_sequence
                    })
            
        except Exception as e:
            raise ValueError(f"Error parsing {file_format} file: {str(e)}")
        
        return sequences
    
    def get_reading_frame_sequence(self, sequence, frame):
        """Get sequence for specific reading frame"""
        if frame.startswith('+'):
            # Forward frames
            frame_num = int(frame[1:]) - 1
            return sequence[frame_num:]
        else:
            # Reverse frames
            frame_num = int(frame[1:]) - 1
            # Get reverse complement
            seq_obj = Seq(sequence)
            reverse_comp = str(seq_obj.reverse_complement())
            return reverse_comp[frame_num:]
    
    def translate_to_mrna(self, dna_sequence):
        """Convert DNA to mRNA (T -> U)"""
        return dna_sequence.replace('T', 'U')
    
    def translate_to_trna(self, mrna_sequence):
        """Convert mRNA to tRNA anticodon sequence"""
        # tRNA anticodon is complementary to mRNA codon
        complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement_map.get(base, base) for base in mrna_sequence)
    
    def translate_to_protein(self, dna_sequence):
        """Translate DNA to protein sequence"""
        seq_obj = Seq(dna_sequence)
        return str(seq_obj.translate())
    
    def translate_sequence(self, sequence, reading_frame, output_types):
        """Main translation function"""
        results = {}
        
        # Get sequence for the specified reading frame
        frame_sequence = self.get_reading_frame_sequence(sequence, reading_frame)
        
        # Ensure sequence length is multiple of 3 for protein translation
        if len(frame_sequence) % 3 != 0:
            frame_sequence = frame_sequence[:-(len(frame_sequence) % 3)]
        
        if 'mrna' in output_types:
            results['mrna'] = self.translate_to_mrna(frame_sequence)
        
        if 'trna' in output_types:
            mrna = self.translate_to_mrna(frame_sequence)
            results['trna'] = self.translate_to_trna(mrna)
        
        if 'protein' in output_types:
            results['protein'] = self.translate_to_protein(frame_sequence)
        
        return results

def main():
    """Main function for command-line execution"""
    parser = argparse.ArgumentParser(description='DNA Translator')
    parser.add_argument('--sequence', required=True, help='DNA sequence to translate')
    parser.add_argument('--reading-frame', default='+1', help='Reading frame (+1, +2, +3, -1, -2, -3)')
    parser.add_argument('--output-types', required=True, help='Comma-separated output types (mrna,trna,protein)')
    parser.add_argument('--file-content', help='File content for file-based translation')
    parser.add_argument('--file-format', default='txt', help='File format (fasta, genbank, txt)')
    
    args = parser.parse_args()
    
    translator = DNATranslatorDirect()
    
    try:
        # Parse output types
        output_types = [t.strip() for t in args.output_types.split(',')]
        
        # Handle file content if provided
        if args.file_content:
            sequences = translator.parse_file_content(args.file_content, args.file_format)
            if not sequences:
                raise ValueError("No valid sequences found in file")
            
            # Use first sequence
            sequence = sequences[0]['sequence']
        else:
            sequence = args.sequence
        
        # Validate DNA sequence
        is_valid, clean_sequence = translator.validate_dna_sequence(sequence)
        if not is_valid:
            raise ValueError("Invalid DNA sequence. Only A, T, G, C, N nucleotides are allowed.")
        
        # Perform translation
        results = translator.translate_sequence(clean_sequence, args.reading_frame, output_types)
        
        # Add sequence info
        seq_obj = Seq(clean_sequence)
        sequence_info = {
            'length': len(clean_sequence),
            'gc_content': round(gc_fraction(seq_obj) * 100, 2),
            'reading_frame': args.reading_frame
        }
        
        # Output JSON result
        output = {
            'success': True,
            'sequence_info': sequence_info,
            'results': results,
            'original_sequence': clean_sequence  # Include original sequence for file parsing
        }
        
        print(json.dumps(output))
        
    except Exception as e:
        # Output error as JSON
        error_output = {
            'success': False,
            'error': str(e)
        }
        print(json.dumps(error_output))
        sys.exit(1)

if __name__ == '__main__':
    main()