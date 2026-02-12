#!/usr/bin/env python3
"""
Visual Codon Finder - Find and highlight specific amino acid codons in DNA sequences
"""

import argparse
import json
import sys
import re
from typing import Dict, List, Tuple, Optional

class CodonFinder:
    def __init__(self):
        """Initialize the Codon Finder with genetic code mapping."""
        # Standard genetic code (DNA codons to amino acids)
        self.genetic_code = {
            # Alanine (A)
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            # Arginine (R)
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
            # Asparagine (N)
            'AAT': 'N', 'AAC': 'N',
            # Aspartic acid (D)
            'GAT': 'D', 'GAC': 'D',
            # Cysteine (C)
            'TGT': 'C', 'TGC': 'C',
            # Glutamic acid (E)
            'GAA': 'E', 'GAG': 'E',
            # Glutamine (Q)
            'CAA': 'Q', 'CAG': 'Q',
            # Glycine (G)
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            # Histidine (H)
            'CAT': 'H', 'CAC': 'H',
            # Isoleucine (I)
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
            # Leucine (L)
            'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            # Lysine (K)
            'AAA': 'K', 'AAG': 'K',
            # Methionine (M) - START
            'ATG': 'M',
            # Phenylalanine (F)
            'TTT': 'F', 'TTC': 'F',
            # Proline (P)
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            # Serine (S)
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
            # Threonine (T)
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            # Tryptophan (W)
            'TGG': 'W',
            # Tyrosine (Y)
            'TAT': 'Y', 'TAC': 'Y',
            # Valine (V)
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            # Stop codons
            'TAA': '*', 'TAG': '*', 'TGA': '*'
        }
        
        # Reverse mapping: amino acid to codons
        self.amino_acid_codons = {}
        for codon, aa in self.genetic_code.items():
            if aa not in self.amino_acid_codons:
                self.amino_acid_codons[aa] = []
            self.amino_acid_codons[aa].append(codon)
        
        # Amino acid information
        self.amino_acid_info = {
            'A': {'name': 'Alanine', 'full_name': 'Alanine', 'symbol': 'Ala', 'type': 'nonpolar'},
            'R': {'name': 'Arginine', 'full_name': 'Arginine', 'symbol': 'Arg', 'type': 'basic'},
            'N': {'name': 'Asparagine', 'full_name': 'Asparagine', 'symbol': 'Asn', 'type': 'polar'},
            'D': {'name': 'Aspartic Acid', 'full_name': 'Aspartic Acid', 'symbol': 'Asp', 'type': 'acidic'},
            'C': {'name': 'Cysteine', 'full_name': 'Cysteine', 'symbol': 'Cys', 'type': 'polar'},
            'E': {'name': 'Glutamic Acid', 'full_name': 'Glutamic Acid', 'symbol': 'Glu', 'type': 'acidic'},
            'Q': {'name': 'Glutamine', 'full_name': 'Glutamine', 'symbol': 'Gln', 'type': 'polar'},
            'G': {'name': 'Glycine', 'full_name': 'Glycine', 'symbol': 'Gly', 'type': 'nonpolar'},
            'H': {'name': 'Histidine', 'full_name': 'Histidine', 'symbol': 'His', 'type': 'basic'},
            'I': {'name': 'Isoleucine', 'full_name': 'Isoleucine', 'symbol': 'Ile', 'type': 'nonpolar'},
            'L': {'name': 'Leucine', 'full_name': 'Leucine', 'symbol': 'Leu', 'type': 'nonpolar'},
            'K': {'name': 'Lysine', 'full_name': 'Lysine', 'symbol': 'Lys', 'type': 'basic'},
            'M': {'name': 'Methionine', 'full_name': 'Methionine', 'symbol': 'Met', 'type': 'nonpolar'},
            'F': {'name': 'Phenylalanine', 'full_name': 'Phenylalanine', 'symbol': 'Phe', 'type': 'nonpolar'},
            'P': {'name': 'Proline', 'full_name': 'Proline', 'symbol': 'Pro', 'type': 'nonpolar'},
            'S': {'name': 'Serine', 'full_name': 'Serine', 'symbol': 'Ser', 'type': 'polar'},
            'T': {'name': 'Threonine', 'full_name': 'Threonine', 'symbol': 'Thr', 'type': 'polar'},
            'W': {'name': 'Tryptophan', 'full_name': 'Tryptophan', 'symbol': 'Trp', 'type': 'nonpolar'},
            'Y': {'name': 'Tyrosine', 'full_name': 'Tyrosine', 'symbol': 'Tyr', 'type': 'polar'},
            'V': {'name': 'Valine', 'full_name': 'Valine', 'symbol': 'Val', 'type': 'nonpolar'},
            '*': {'name': 'Stop', 'full_name': 'Stop Codon', 'symbol': 'Stop', 'type': 'stop'}
        }
        
        # Color mapping for amino acids
        self.amino_acid_colors = {
            'A': '#3B82F6',    # Blue - Alanine
            'R': '#06B6D4',    # Cyan - Arginine  
            'N': '#8B5CF6',    # Purple - Asparagine
            'D': '#EAB308',    # Yellow - Aspartic Acid
            'C': '#EF4444',    # Red - Cysteine
            'E': '#10B981',    # Green - Glutamic Acid
            'Q': '#F97316',    # Orange - Glutamine
            'G': '#EC4899',    # Pink - Glycine
            'H': '#6366F1',    # Indigo - Histidine
            'I': '#14B8A6',    # Teal - Isoleucine
            'L': '#84CC16',    # Lime - Leucine
            'K': '#0EA5E9',    # Sky - Lysine
            'M': '#10B981',    # Emerald - Methionine (START)
            'F': '#A855F7',    # Violet - Phenylalanine
            'P': '#F43F5E',    # Rose - Proline
            'S': '#F59E0B',    # Amber - Serine
            'T': '#059669',    # Emerald - Threonine
            'W': '#D946EF',    # Fuchsia - Tryptophan
            'Y': '#64748B',    # Slate - Tyrosine
            'V': '#71717A',    # Zinc - Valine
            '*': '#DC2626'     # Red - Stop codons
        }
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and validate DNA sequence."""
        # Remove whitespace and convert to uppercase
        cleaned = re.sub(r'\s+', '', sequence.upper())
        
        # Remove non-DNA characters
        cleaned = re.sub(r'[^ATCG]', '', cleaned)
        
        return cleaned
    
    def find_codons(self, sequence: str, target_amino_acid: str) -> Dict:
        """Find all codons for a specific amino acid in the sequence."""
        
        try:
            # Clean sequence
            clean_seq = self.clean_sequence(sequence)
            
            if len(clean_seq) < 3:
                return {
                    'success': False,
                    'error': 'Sequence too short (minimum 3 nucleotides required)'
                }
            
            # Get target amino acid info
            if target_amino_acid not in self.amino_acid_info:
                return {
                    'success': False,
                    'error': f'Unknown amino acid: {target_amino_acid}'
                }
            
            # Get codons for target amino acid
            target_codons = self.amino_acid_codons.get(target_amino_acid, [])
            
            # Find all codons in sequence
            all_codons = []
            target_positions = []
            
            # Process sequence in triplets
            for i in range(0, len(clean_seq) - 2, 3):
                codon = clean_seq[i:i+3]
                amino_acid = self.genetic_code.get(codon, 'X')  # X for unknown
                
                codon_info = {
                    'position': i,
                    'codon': codon,
                    'amino_acid': amino_acid,
                    'is_target': amino_acid == target_amino_acid,
                    'color': self.amino_acid_colors.get(amino_acid, '#6B7280')
                }
                
                all_codons.append(codon_info)
                
                if amino_acid == target_amino_acid:
                    target_positions.append(i)
            
            # Create highlighted sequence
            highlighted_sequence = self.create_highlighted_sequence(clean_seq, all_codons, target_amino_acid)
            
            # Calculate statistics
            total_codons = len(all_codons)
            target_count = len(target_positions)
            frequency = (target_count / total_codons * 100) if total_codons > 0 else 0
            
            # Get amino acid details
            aa_info = self.amino_acid_info[target_amino_acid]
            
            return {
                'success': True,
                'sequence': clean_seq,
                'target_amino_acid': target_amino_acid,
                'amino_acid_info': aa_info,
                'target_codons': target_codons,
                'all_codons': all_codons,
                'target_positions': target_positions,
                'highlighted_sequence': highlighted_sequence,
                'statistics': {
                    'sequence_length': len(clean_seq),
                    'total_codons': total_codons,
                    'target_count': target_count,
                    'frequency': round(frequency, 2)
                },
                'codon_details': self.generate_codon_details(target_positions, clean_seq, target_amino_acid)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f'Analysis failed: {str(e)}'
            }
    
    def create_highlighted_sequence(self, sequence: str, all_codons: List[Dict], target_amino_acid: str) -> str:
        """Create HTML highlighted sequence."""
        highlighted = ""
        
        for i, codon_info in enumerate(all_codons):
            codon = codon_info['codon']
            amino_acid = codon_info['amino_acid']
            is_target = codon_info['is_target']
            color = codon_info['color']
            
            # Create tooltip info
            aa_info = self.amino_acid_info.get(amino_acid, {})
            tooltip_text = f"{aa_info.get('full_name', 'Unknown')} ({aa_info.get('symbol', 'X')}) - {codon}"
            
            if is_target:
                # Highlight target codons
                highlighted += f'<span class="codon-highlight target-codon" style="background-color: {color}20; color: {color}; border: 1px solid {color}50; padding: 1px 2px; border-radius: 3px; font-weight: bold;" title="{tooltip_text}" data-codon="{codon}" data-amino="{amino_acid}" data-position="{codon_info["position"]}">{codon}</span>'
            else:
                # Regular codons
                highlighted += f'<span class="codon-normal" style="color: #9CA3AF;" title="{tooltip_text}" data-codon="{codon}" data-amino="{amino_acid}" data-position="{codon_info["position"]}">{codon}</span>'
            
            # Add space between codons
            if i < len(all_codons) - 1:
                highlighted += " "
        
        return highlighted
    
    def generate_codon_details(self, target_positions: List[int], sequence: str, target_amino_acid: str) -> List[Dict]:
        """Generate detailed information about found codons."""
        details = []
        
        aa_info = self.amino_acid_info[target_amino_acid]
        
        for i, pos in enumerate(target_positions):
            codon = sequence[pos:pos+3]
            
            details.append({
                'index': i + 1,
                'position': pos + 1,  # 1-based position
                'codon': codon,
                'amino_acid': target_amino_acid,
                'amino_acid_name': aa_info['full_name'],
                'amino_acid_symbol': aa_info['symbol'],
                'amino_acid_type': aa_info['type']
            })
        
        return details

def main():
    parser = argparse.ArgumentParser(description='Visual Codon Finder - Find and highlight amino acid codons')
    parser.add_argument('--sequence', help='DNA sequence to analyze')
    parser.add_argument('--amino-acid', help='Target amino acid (single letter code)')
    
    args = parser.parse_args()
    
    # Get sequence and amino acid from arguments or stdin
    sequence = args.sequence
    amino_acid = args.amino_acid
    
    # If no sequence from args, try stdin
    if not sequence or not amino_acid:
        try:
            input_data = sys.stdin.read()
            if input_data:
                data = json.loads(input_data)
                sequence = sequence or data.get('sequence', '')
                amino_acid = amino_acid or data.get('amino_acid', '')
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
    
    if not sequence or not amino_acid:
        print(json.dumps({
            'success': False,
            'error': 'Sequence and amino acid are required'
        }))
        return
    
    finder = CodonFinder()
    
    # Find codons
    result = finder.find_codons(sequence.upper(), amino_acid.upper())
    
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()