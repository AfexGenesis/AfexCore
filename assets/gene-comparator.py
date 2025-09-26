#!/usr/bin/env python3
"""
AFEX Gene Comparator - DIAMOND-Inspired Tool
User-friendly pairwise gene comparison with detailed statistics
"""

import sys
import json
import argparse
import time
from io import StringIO
import math

try:
    from Bio import pairwise2
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction
    from Bio import SeqIO
    from Bio.pairwise2 import format_alignment
    from Bio.Align import PairwiseAligner
except ImportError as e:
    print(json.dumps({
        'success': False,
        'error': f'Missing required library: {str(e)}. Please install: pip install biopython'
    }))
    sys.exit(1)

class AFEXGeneComparator:
    """DIAMOND-inspired gene comparison tool"""
    
    def __init__(self):
        # Scoring parameters (DIAMOND-like)
        self.match_score = 2
        self.mismatch_score = -1
        self.gap_open = -2
        self.gap_extend = -0.5
        
        # Statistical parameters
        self.lambda_param = 1.33
        self.k_param = 0.621
        
    def parse_sequence_input(self, sequence_data, sequence_name="sequence"):
        """Parse sequence from various formats (FASTA, GenBank, or plain text)"""
        sequences = []
        
        try:
            # Try to parse as FASTA first
            if sequence_data.strip().startswith('>'):
                fasta_io = StringIO(sequence_data)
                for record in SeqIO.parse(fasta_io, "fasta"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq).upper()
                    })
            
            # Try GenBank format
            elif 'LOCUS' in sequence_data and 'ORIGIN' in sequence_data:
                genbank_io = StringIO(sequence_data)
                for record in SeqIO.parse(genbank_io, "genbank"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq).upper()
                    })
            
            # Plain text sequence
            else:
                clean_sequence = ''.join(c.upper() for c in sequence_data if c.upper() in 'ATGCN')
                if clean_sequence:
                    sequences.append({
                        'id': sequence_name,
                        'description': f'{sequence_name} (plain text)',
                        'sequence': clean_sequence
                    })
            
            return sequences
            
        except Exception as e:
            raise ValueError(f"Error parsing sequence: {str(e)}")
    
    def validate_dna_sequence(self, sequence):
        """Validate DNA sequence contains only valid nucleotides"""
        valid_bases = set('ATGCN')
        sequence_upper = sequence.upper()
        invalid_bases = set(sequence_upper) - valid_bases
        
        if invalid_bases:
            return False, f"Invalid bases found: {', '.join(invalid_bases)}"
        
        return True, sequence_upper
    
    def perform_alignment(self, seq1, seq2):
        """Perform pairwise alignment using BioPython"""
        try:
            # Use PairwiseAligner for better control
            aligner = PairwiseAligner()
            aligner.match_score = self.match_score
            aligner.mismatch_score = self.mismatch_score
            aligner.open_gap_score = self.gap_open
            aligner.extend_gap_score = self.gap_extend
            aligner.mode = 'local'  # Local alignment like DIAMOND
            
            # Perform alignment
            alignments = aligner.align(seq1, seq2)
            best_alignment = alignments[0]  # Get best alignment
            
            return best_alignment
            
        except Exception as e:
            # Fallback to pairwise2 if PairwiseAligner fails
            alignments = pairwise2.align.localms(
                seq1, seq2,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend,
                one_alignment_only=True
            )
            
            if alignments:
                return alignments[0]
            else:
                raise ValueError("Alignment failed")
    
    def calculate_statistics(self, alignment, seq1, seq2):
        """Calculate DIAMOND-like statistics from alignment"""
        try:
            # Handle different alignment object types
            if hasattr(alignment, 'aligned'):
                # New PairwiseAligner format
                aligned_seq1 = str(alignment.query)
                aligned_seq2 = str(alignment.target)
                score = alignment.score
            else:
                # Old pairwise2 format
                aligned_seq1 = alignment[0]
                aligned_seq2 = alignment[1]
                score = alignment[2]
            
            # Calculate basic statistics
            alignment_length = len(aligned_seq1)
            matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
            mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
            gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
            
            # Count gap opens (consecutive gaps count as one open)
            gap_opens = 0
            in_gap = False
            for a, b in zip(aligned_seq1, aligned_seq2):
                if a == '-' or b == '-':
                    if not in_gap:
                        gap_opens += 1
                        in_gap = True
                else:
                    in_gap = False
            
            # Calculate percentages
            identity_percent = (matches / alignment_length * 100) if alignment_length > 0 else 0
            
            # Calculate similarity (matches + conservative substitutions)
            # For DNA, we'll consider purines (A,G) and pyrimidines (C,T) as similar
            similar_pairs = 0
            for a, b in zip(aligned_seq1, aligned_seq2):
                if a == b and a != '-':
                    similar_pairs += 1
                elif (a in 'AG' and b in 'AG') or (a in 'CT' and b in 'CT'):
                    similar_pairs += 1
            
            similarity_percent = (similar_pairs / alignment_length * 100) if alignment_length > 0 else 0
            
            # Calculate coverage
            query_coverage = (alignment_length - aligned_seq1.count('-')) / len(seq1) * 100
            subject_coverage = (alignment_length - aligned_seq2.count('-')) / len(seq2) * 100
            
            # Calculate E-value (simplified approximation)
            effective_length = len(seq1) * len(seq2)
            bit_score = (self.lambda_param * score - math.log(self.k_param)) / math.log(2)
            e_value = effective_length * math.pow(2, -bit_score) if bit_score > 0 else 1.0
            
            return {
                'alignment_length': alignment_length,
                'matches': matches,
                'mismatches': mismatches,
                'gaps': gaps,
                'gap_opens': gap_opens,
                'identity_percent': round(identity_percent, 2),
                'similarity_percent': round(similarity_percent, 2),
                'query_coverage': round(query_coverage, 2),
                'subject_coverage': round(subject_coverage, 2),
                'score': round(score, 1),
                'bit_score': round(bit_score, 1),
                'e_value': f"{e_value:.2e}",
                'aligned_query': aligned_seq1,
                'aligned_subject': aligned_seq2
            }
            
        except Exception as e:
            raise ValueError(f"Error calculating statistics: {str(e)}")
    
    def analyze_sequences(self, seq1_info, seq2_info):
        """Analyze sequence composition and properties"""
        seq1 = seq1_info['sequence']
        seq2 = seq2_info['sequence']
        
        # Basic composition
        def get_composition(seq):
            total = len(seq)
            return {
                'length': total,
                'gc_content': round(gc_fraction(Seq(seq)) * 100, 2),
                'a_count': seq.count('A'),
                't_count': seq.count('T'),
                'g_count': seq.count('G'),
                'c_count': seq.count('C'),
                'n_count': seq.count('N')
            }
        
        return {
            'query': get_composition(seq1),
            'subject': get_composition(seq2)
        }
    
    def format_alignment_display(self, aligned_seq1, aligned_seq2, line_length=60):
        """Format alignment for display with match indicators"""
        lines = []
        
        for i in range(0, len(aligned_seq1), line_length):
            query_line = aligned_seq1[i:i+line_length]
            subject_line = aligned_seq2[i:i+line_length]
            
            # Create match line
            match_line = ''
            for q, s in zip(query_line, subject_line):
                if q == s and q != '-':
                    match_line += '|'
                elif q != '-' and s != '-':
                    match_line += ':'
                else:
                    match_line += ' '
            
            lines.append(f"Query  {i+1:>6} {query_line} {i+len(query_line)}")
            lines.append(f"              {match_line}")
            lines.append(f"Sbjct  {i+1:>6} {subject_line} {i+len(subject_line)}")
            lines.append("")
        
        return '\n'.join(lines)
    
    def compare_genes(self, seq1_data, seq2_data):
        """Main comparison function"""
        try:
            # Parse sequences
            seq1_list = self.parse_sequence_input(seq1_data, "Query")
            seq2_list = self.parse_sequence_input(seq2_data, "Subject")
            
            if not seq1_list or not seq2_list:
                raise ValueError("Could not parse input sequences")
            
            # Use first sequence from each input
            seq1_info = seq1_list[0]
            seq2_info = seq2_list[0]
            
            # Validate sequences
            is_valid1, seq1_clean = self.validate_dna_sequence(seq1_info['sequence'])
            is_valid2, seq2_clean = self.validate_dna_sequence(seq2_info['sequence'])
            
            if not is_valid1:
                raise ValueError(f"Invalid query sequence: {seq1_clean}")
            if not is_valid2:
                raise ValueError(f"Invalid subject sequence: {seq2_clean}")
            
            # Update with clean sequences
            seq1_info['sequence'] = seq1_clean
            seq2_info['sequence'] = seq2_clean
            
            # Perform alignment
            alignment = self.perform_alignment(seq1_clean, seq2_clean)
            
            # Calculate statistics
            stats = self.calculate_statistics(alignment, seq1_clean, seq2_clean)
            
            # Analyze sequence composition
            composition = self.analyze_sequences(seq1_info, seq2_info)
            
            # Format alignment display
            alignment_display = self.format_alignment_display(
                stats['aligned_query'], 
                stats['aligned_subject']
            )
            
            return {
                'success': True,
                'query_info': seq1_info,
                'subject_info': seq2_info,
                'statistics': stats,
                'composition': composition,
                'alignment_display': alignment_display,
                'summary': self.generate_summary(stats, seq1_info, seq2_info)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }
    
    def generate_summary(self, stats, seq1_info, seq2_info):
        """Generate a human-readable summary"""
        identity = stats['identity_percent']
        
        if identity >= 95:
            similarity_desc = "Nearly identical genes"
        elif identity >= 80:
            similarity_desc = "Highly similar genes"
        elif identity >= 60:
            similarity_desc = "Moderately similar genes"
        elif identity >= 40:
            similarity_desc = "Distantly related genes"
        else:
            similarity_desc = "Low similarity genes"
        
        return {
            'description': similarity_desc,
            'identity_category': 'high' if identity >= 80 else 'medium' if identity >= 60 else 'low',
            'functional_prediction': self.predict_function(stats)
        }
    
    def predict_function(self, stats):
        """Predict functional relationship based on similarity"""
        identity = stats['identity_percent']
        coverage = min(stats['query_coverage'], stats['subject_coverage'])
        
        if identity >= 95 and coverage >= 90:
            return "Likely identical genes or very close orthologs"
        elif identity >= 80 and coverage >= 70:
            return "Probable orthologs with conserved function"
        elif identity >= 60 and coverage >= 50:
            return "Possible paralogs or distant orthologs"
        elif identity >= 40:
            return "Distantly related, may share common ancestry"
        else:
            return "Low similarity, unlikely to share function"

def main():
    """Main function for command-line execution"""
    parser = argparse.ArgumentParser(description='AFEX Gene Comparator - DIAMOND-inspired tool')
    parser.add_argument('--query', required=True, help='First gene sequence (query)')
    parser.add_argument('--subject', required=True, help='Second gene sequence (subject)')
    parser.add_argument('--format', default='auto', help='Input format (auto, fasta, genbank, text)')
    
    args = parser.parse_args()
    
    comparator = AFEXGeneComparator()
    
    try:
        # Check if we need to read from stdin
        if args.query == 'STDIN_QUERY' and args.subject == 'STDIN_SUBJECT':
            # Read JSON input from stdin
            stdin_data = sys.stdin.read()
            try:
                input_json = json.loads(stdin_data)
                query_sequence = input_json['query']
                subject_sequence = input_json['subject']
            except (json.JSONDecodeError, KeyError) as e:
                raise ValueError(f"Invalid JSON input from stdin: {str(e)}")
        else:
            # Use command line arguments
            query_sequence = args.query
            subject_sequence = args.subject
        
        # Perform comparison
        result = comparator.compare_genes(query_sequence, subject_sequence)
        
        # Output JSON result
        print(json.dumps(result, indent=2))
        
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