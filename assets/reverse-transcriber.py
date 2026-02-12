#!/usr/bin/env python3
"""
AFEX Genesis™ - Reverse Transcriber (No API)
Command-line reverse transcription: mRNA/tRNA/Protein → DNA
"""

import sys
import json
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import io
import re

class ReverseTranscriber:
    """Reverse transcription from mRNA/tRNA/Protein to DNA"""
    
    def __init__(self):
        # Standard genetic code (codon to amino acid)
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
        
        # Reverse genetic code (amino acid to possible codons)
        self.reverse_genetic_code = {}
        for codon, aa in self.genetic_code.items():
            if aa not in self.reverse_genetic_code:
                self.reverse_genetic_code[aa] = []
            self.reverse_genetic_code[aa].append(codon)
        
        # Most common codons for each amino acid (for optimization)
        self.preferred_codons = {
            'F': 'TTT', 'L': 'CTG', 'S': 'TCT', 'Y': 'TAT',
            'C': 'TGT', 'W': 'TGG', 'P': 'CCT', 'H': 'CAT',
            'Q': 'CAG', 'R': 'CGT', 'I': 'ATT', 'M': 'ATG',
            'T': 'ACT', 'N': 'AAT', 'K': 'AAG', 'V': 'GTT',
            'A': 'GCT', 'D': 'GAT', 'E': 'GAG', 'G': 'GGT',
            '*': 'TAA'  # Stop codon
        }
    
    def validate_sequence(self, sequence, seq_type):
        """Validate input sequence based on type"""
        sequence = sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '')
        
        if seq_type == 'mrna':
            valid_bases = set('AUGC')
            if not all(base in valid_bases for base in sequence):
                raise ValueError("mRNA sequence contains invalid bases. Only A, U, G, C allowed.")
        
        elif seq_type == 'trna':
            valid_bases = set('AUGC')
            if not all(base in valid_bases for base in sequence):
                raise ValueError("tRNA sequence contains invalid bases. Only A, U, G, C allowed.")
        
        elif seq_type == 'protein':
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
            if not all(aa in valid_aa for aa in sequence):
                raise ValueError("Protein sequence contains invalid amino acids.")
        
        return sequence
    
    def mrna_to_dna(self, mrna_sequence):
        """Convert mRNA to DNA"""
        # mRNA is transcribed from DNA template strand
        # So mRNA sequence is same as DNA coding strand (except U→T)
        dna_coding = mrna_sequence.replace('U', 'T')
        
        # Template strand is complement of coding strand
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        dna_template = ''.join(complement_map[base] for base in dna_coding)
        
        return {
            'template_strand': dna_template,
            'coding_strand': dna_coding,
            'primary_result': dna_coding  # Return coding strand as primary
        }
    
    def trna_to_dna(self, trna_sequence):
        """Convert tRNA to DNA"""
        # tRNA is complementary to mRNA, so first convert to mRNA
        # tRNA → mRNA (complement)
        rna_complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        mrna_sequence = ''.join(rna_complement_map[base] for base in trna_sequence)
        
        # Then mRNA → DNA
        return self.mrna_to_dna(mrna_sequence)
    
    def protein_to_dna(self, protein_sequence, codon_usage='optimal'):
        """Convert protein sequence to DNA using reverse translation"""
        dna_sequence = ""
        
        for amino_acid in protein_sequence:
            if amino_acid not in self.reverse_genetic_code:
                raise ValueError(f"Unknown amino acid: {amino_acid}")
            
            # Choose codon based on usage preference
            if codon_usage == 'optimal' and amino_acid in self.preferred_codons:
                codon = self.preferred_codons[amino_acid]
            else:
                # Use first available codon
                codon = self.reverse_genetic_code[amino_acid][0]
            
            dna_sequence += codon
        
        return {
            'coding_strand': dna_sequence,
            'template_strand': self.get_complement(dna_sequence)[::-1],
            'primary_result': dna_sequence,
            'alternative_sequences': self.get_alternative_sequences(protein_sequence)
        }
    
    def get_complement(self, dna_sequence):
        """Get DNA complement"""
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement_map[base] for base in dna_sequence)
    
    def get_alternative_sequences(self, protein_sequence, max_alternatives=3):
        """Generate alternative DNA sequences for the same protein"""
        alternatives = []
        
        # Generate a few alternative sequences using different codon choices
        for alt_num in range(min(max_alternatives, 3)):
            alt_sequence = ""
            for i, amino_acid in enumerate(protein_sequence):
                available_codons = self.reverse_genetic_code[amino_acid]
                # Use different codon for variety
                codon_index = (alt_num + i) % len(available_codons)
                alt_sequence += available_codons[codon_index]
            alternatives.append(alt_sequence)
        
        return alternatives
    
    def parse_file_content(self, file_content, file_format):
        """Parse file content and extract sequence"""
        sequences = []
        
        if file_format.lower() == 'fasta':
            # Parse FASTA format
            fasta_io = io.StringIO(file_content)
            for record in SeqIO.parse(fasta_io, "fasta"):
                sequences.append(str(record.seq))
        
        elif file_format.lower() == 'txt':
            # Simple text format - extract sequence lines
            lines = file_content.strip().split('\n')
            sequence = ""
            for line in lines:
                line = line.strip()
                if line and not line.startswith('>') and not line.startswith('#'):
                    sequence += line
            if sequence:
                sequences.append(sequence)
        
        else:
            # Default: treat as plain sequence
            clean_sequence = re.sub(r'[^A-Za-z*]', '', file_content)
            if clean_sequence:
                sequences.append(clean_sequence)
        
        return sequences[0] if sequences else ""
    
    def reverse_transcribe(self, sequence, input_type, codon_usage='optimal'):
        """Main reverse transcription function"""
        
        # Validate and clean sequence
        clean_sequence = self.validate_sequence(sequence, input_type)
        
        if not clean_sequence:
            raise ValueError("Empty sequence provided")
        
        # Perform reverse transcription based on input type
        if input_type == 'mrna':
            result = self.mrna_to_dna(clean_sequence)
            result['input_type'] = 'mRNA'
            result['process'] = 'mRNA → DNA'
            
        elif input_type == 'trna':
            result = self.trna_to_dna(clean_sequence)
            result['input_type'] = 'tRNA'
            result['process'] = 'tRNA → mRNA → DNA'
            
        elif input_type == 'protein':
            result = self.protein_to_dna(clean_sequence, codon_usage)
            result['input_type'] = 'Protein'
            result['process'] = 'Protein → DNA (reverse translation)'
            
        else:
            raise ValueError(f"Unsupported input type: {input_type}")
        
        return result

def main():
    parser = argparse.ArgumentParser(description='AFEX Genesis™ Reverse Transcriber')
    parser.add_argument('--sequence', help='Input sequence')
    parser.add_argument('--input-type', choices=['mrna', 'trna', 'protein'], 
                       help='Type of input sequence')
    parser.add_argument('--codon-usage', default='optimal', choices=['optimal', 'first'],
                       help='Codon usage preference for protein reverse translation')
    parser.add_argument('--output-types', default='dna', 
                       help='Comma-separated output types (dna,template,alternatives)')
    parser.add_argument('--file-content', help='File content for parsing')
    parser.add_argument('--file-format', help='File format (fasta, txt)')
    
    args = parser.parse_args()
    
    # Get inputs from arguments or stdin
    sequence = args.sequence
    input_type = args.input_type
    codon_usage = args.codon_usage
    output_types_arg = args.output_types
    file_content = args.file_content
    file_format = args.file_format
    
    # If no sequence from args, try stdin
    if not sequence or not input_type:
        try:
            input_data = sys.stdin.read()
            if input_data:
                data = json.loads(input_data)
                sequence = sequence or data.get('sequence', '')
                input_type = input_type or data.get('input_type', '')
                codon_usage = data.get('codon_usage', codon_usage)
                output_types_arg = data.get('output_types', output_types_arg)
                file_content = data.get('file_content', file_content)
                file_format = data.get('file_format', file_format)
        except json.JSONDecodeError:
            error_output = {
                'success': False,
                'error': 'Failed to parse JSON from stdin'
            }
            print(json.dumps(error_output))
            sys.exit(1)
        except Exception as e:
            error_output = {
                'success': False,
                'error': f'Error reading from stdin: {str(e)}'
            }
            print(json.dumps(error_output))
            sys.exit(1)
    
    if not sequence or not input_type:
        error_output = {
            'success': False,
            'error': 'Sequence and input type are required'
        }
        print(json.dumps(error_output))
        sys.exit(1)
    
    try:
        transcriber = ReverseTranscriber()
        
        # Handle file input if provided
        if file_content and file_format:
            sequence = transcriber.parse_file_content(file_content, file_format)
            if not sequence:
                raise ValueError("No valid sequence found in file")
        
        # Perform reverse transcription
        result = transcriber.reverse_transcribe(sequence, input_type, codon_usage)
        
        # Prepare output based on requested types
        output_types = [t.strip() for t in output_types_arg.split(',')]
        results = {}
        
        if 'dna' in output_types:
            results['dna'] = result['primary_result']
        
        if 'template' in output_types and 'template_strand' in result:
            results['template'] = result['template_strand']
        
        if 'alternatives' in output_types and 'alternative_sequences' in result:
            results['alternatives'] = result['alternative_sequences'][:3]  # Limit to 3
        
        # Add sequence info
        seq_obj = Seq(result['primary_result'])
        sequence_info = {
            'input_length': len(sequence),
            'output_length': len(result['primary_result']),
            'gc_content': round(gc_fraction(seq_obj) * 100, 2),
            'input_type': result['input_type'],
            'process': result['process']
        }
        
        # Output JSON result
        output = {
            'success': True,
            'sequence_info': sequence_info,
            'results': results,
            'original_sequence': sequence,
            'primary_dna': result['primary_result']
        }
        
        print(json.dumps(output))
        
    except Exception as e:
        error_output = {
            'success': False,
            'error': str(e)
        }
        print(json.dumps(error_output))
        sys.exit(1)

if __name__ == "__main__":
    main()