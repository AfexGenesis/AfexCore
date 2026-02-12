#!/usr/bin/env python3
"""
Base Counter - Simple DNA sequence analysis tool
Counts nucleotide bases and calculates GC content
"""

import argparse
import json
import sys
import re
from typing import Dict

class BaseCounter:
    def __init__(self):
        """Initialize the Base Counter."""
        pass
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and validate DNA sequence."""
        # Remove whitespace and convert to uppercase
        sequence = re.sub(r'\s+', '', sequence.upper())
        
        # Remove non-DNA characters
        sequence = re.sub(r'[^ATGC]', '', sequence)
        
        return sequence
    
    def count_bases(self, sequence: str) -> Dict:
        """Count individual bases in the sequence."""
        clean_seq = self.clean_sequence(sequence)
        
        if not clean_seq:
            return {
                'success': False,
                'error': 'Invalid or empty DNA sequence'
            }
        
        # Count each base
        base_counts = {
            'A': clean_seq.count('A'),
            'T': clean_seq.count('T'),
            'G': clean_seq.count('G'),
            'C': clean_seq.count('C')
        }
        
        total_bases = len(clean_seq)
        
        # Calculate percentages
        base_percentages = {}
        for base, count in base_counts.items():
            percentage = (count / total_bases * 100) if total_bases > 0 else 0
            base_percentages[base] = round(percentage, 2)
        
        # Calculate GC content
        gc_count = base_counts['G'] + base_counts['C']
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
        
        # Calculate AT content
        at_count = base_counts['A'] + base_counts['T']
        at_content = (at_count / total_bases * 100) if total_bases > 0 else 0
        
        # Classify GC content
        gc_classification = self.classify_gc_content(gc_content)
        
        # Estimate melting temperature (simple approximation)
        melting_temp = self.estimate_melting_temperature(clean_seq)
        
        return {
            'success': True,
            'original_sequence': sequence,
            'clean_sequence': clean_seq,
            'total_bases': total_bases,
            'base_counts': base_counts,
            'base_percentages': base_percentages,
            'gc_content': round(gc_content, 2),
            'at_content': round(at_content, 2),
            'gc_classification': gc_classification,
            'melting_temperature': round(melting_temp, 1)
        }
    
    def classify_gc_content(self, gc_content: float) -> str:
        """Classify GC content into categories."""
        if gc_content < 30:
            return "Low GC (AT-rich)"
        elif gc_content < 45:
            return "Moderate Low GC"
        elif gc_content < 55:
            return "Balanced GC"
        elif gc_content < 65:
            return "Moderate High GC"
        else:
            return "High GC (GC-rich)"
    
    def estimate_melting_temperature(self, sequence: str) -> float:
        """Estimate melting temperature using simple approximation."""
        if len(sequence) == 0:
            return 0.0
        
        # Simple approximation: Tm = 2(A+T) + 4(G+C)
        # This is a very basic formula, more accurate for short sequences
        a_count = sequence.count('A')
        t_count = sequence.count('T')
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        
        if len(sequence) <= 14:
            # For short sequences (primers)
            tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
        else:
            # For longer sequences, use a different approximation
            gc_content = ((g_count + c_count) / len(sequence)) * 100
            tm = 81.5 + 0.41 * gc_content - (675 / len(sequence))
        
        return max(0, tm)  # Ensure non-negative temperature

def main():
    parser = argparse.ArgumentParser(description='Base Counter - Analyze DNA sequence composition')
    parser.add_argument('--sequence', help='DNA sequence to analyze')
    parser.add_argument('--file-content', help='File content for parsing')
    parser.add_argument('--file-format', help='File format (fasta or txt)')
    
    args = parser.parse_args()
    
    counter = BaseCounter()
    
    # Handle input from stdin if no sequence argument provided
    sequence = args.sequence
    if not sequence:
        try:
            input_data = sys.stdin.read()
            if input_data:
                data = json.loads(input_data)
                sequence = data.get('sequence', '')
        except json.JSONDecodeError:
            print(json.dumps({
                'success': False,
                'error': 'Failed to parse JSON from stdin'
            }))
            return
        except Exception as e:
            print(json.dumps({
                'success': False,
                'error': f'Error reading from stdin: {str(e)}'
            }))
            return
    
    # Handle file input
    if args.file_content and args.file_format:
        try:
            if args.file_format.lower() == 'fasta':
                # Parse FASTA format
                lines = args.file_content.strip().split('\n')
                sequence_lines = [line for line in lines if not line.startswith('>')]
                sequence = ''.join(sequence_lines)
            else:
                # Plain text
                sequence = args.file_content.strip()
        except Exception as e:
            print(json.dumps({
                'success': False,
                'error': f'Failed to parse file: {str(e)}'
            }))
            return
    
    # Analyze sequence
    result = counter.count_bases(sequence)
    
    print(json.dumps(result))

if __name__ == '__main__':
    main()