#!/usr/bin/env python3
"""
AFEX Genesis™ - DNA Sequencer
Advanced DNA sequence analysis and visualization
"""

import sys
import json
import argparse
import re
from io import StringIO
from typing import Dict, List, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction
    from Bio.SeqRecord import SeqRecord
except ImportError as e:
    print(json.dumps({
        'success': False,
        'error': f'Missing required library: {str(e)}. Please install: pip install biopython'
    }))
    sys.exit(1)


class DNASequencer:
    """Advanced DNA Sequencer with comprehensive analysis"""
    
    def __init__(self):
        # Comprehensive restriction enzyme database
        self.restriction_enzymes = {
            # Common 6-base cutters
            'BamHI': 'GGATCC', 'EcoRI': 'GAATTC', 'HindIII': 'AAGCTT', 'XbaI': 'TCTAGA',
            'SalI': 'GTCGAC', 'PstI': 'CTGCAG', 'SmaI': 'CCCGGG', 'KpnI': 'GGTACC',
            'SacI': 'GAGCTC', 'XhoI': 'CTCGAG', 'NotI': 'GCGGCCGC', 'ApaI': 'GGGCCC',
            'BglII': 'AGATCT', 'ClaI': 'ATCGAT', 'NcoI': 'CCATGG', 'NdeI': 'CATATG',
            'NheI': 'GCTAGC', 'SpeI': 'ACTAGT', 'AvrII': 'CCTAGG', 'NsiI': 'ATGCAT',
            'EcoRV': 'GATATC', 'PvuII': 'CAGCTG', 'ScaI': 'AGTACT', 'StuI': 'AGGCCT',
            
            # 4-base cutters (frequent)
            'AluI': 'AGCT', 'HaeIII': 'GGCC', 'HhaI': 'GCGC', 'HpaII': 'CCGG',
            'MspI': 'CCGG', 'RsaI': 'GTAC', 'TaqI': 'TCGA', 'MseI': 'TTAA',
            
            # 8-base rare cutters
            'AsiSI': 'GCGATCGC', 'FseI': 'GGCCGGCC', 'PacI': 'TTAATTAA',
            'SbfI': 'CCTGCAGG', 'SwaI': 'ATTTAAAT', 'PmeI': 'GTTTAAAC',
            
            # Blunt cutters
            'HpaI': 'GTTAAC', 'NruI': 'TCGCGA', 'SspI': 'AATATT',
            
            # Methylation sensitive
            'BstUI': 'CGCG', 'AciI': 'CCGC', 'HinP1I': 'GCGC',
            
            # Type IIS enzymes
            'BsaI': 'GGTCTC', 'BsmBI': 'CGTCTC', 'BbsI': 'GAAGAC', 'Esp3I': 'CGTCTC',
            'SapI': 'GCTCTTC', 'BtgZI': 'GCGATG',
            
            # Nicking enzymes
            'Nb.BsrDI': 'GCAATG', 'Nb.BsmI': 'GAATGC', 'Nt.BspQI': 'GCTCTTC', 'Nb.BtsI': 'GCAGTG'
        }
        
        # Color schemes for nucleotide display
        self.color_schemes = {
            'standard': {'A': 'text-red-400', 'T': 'text-blue-400', 'G': 'text-green-400', 'C': 'text-yellow-400'},
            'purine-pyrimidine': {'A': 'text-blue-400', 'G': 'text-blue-600', 'T': 'text-red-400', 'C': 'text-red-600'},
            'gc-at': {'G': 'text-green-400', 'C': 'text-green-600', 'A': 'text-orange-400', 'T': 'text-orange-600'},
            'rainbow': {'A': 'text-red-400', 'T': 'text-orange-400', 'G': 'text-green-400', 'C': 'text-purple-400'},
            'monochrome': {'A': 'text-gray-400', 'T': 'text-gray-500', 'G': 'text-gray-600', 'C': 'text-gray-700'}
        }
    
    def parse_file_content(self, file_content: str, file_format: str) -> List[Dict]:
        """Parse different file formats and extract DNA sequences"""
        sequences = []
        
        try:
            if file_format.lower() in ['fasta', 'fa']:
                # Parse FASTA format
                fasta_io = StringIO(file_content)
                for record in SeqIO.parse(fasta_io, "fasta"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq).upper()
                    })
            
            elif file_format.lower() in ['genbank', 'gb', 'gbk']:
                # Parse GenBank format
                genbank_io = StringIO(file_content)
                for record in SeqIO.parse(genbank_io, "genbank"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq).upper()
                    })
            
            elif file_format.lower() in ['txt', 'seq']:
                # Parse plain text - assume it's just DNA sequence
                clean_sequence = self.clean_sequence(file_content)
                if clean_sequence:
                    sequences.append({
                        'id': 'sequence_1',
                        'description': 'Plain text sequence',
                        'sequence': clean_sequence
                    })
            
        except Exception as e:
            raise ValueError(f"Error parsing {file_format} file: {str(e)}")
        
        return sequences
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and validate DNA sequence"""
        # Remove whitespace, numbers, and convert to uppercase
        sequence = re.sub(r'[\s\d]', '', sequence.upper())
        
        # Remove non-DNA characters (keep only A, T, G, C, N)
        sequence = re.sub(r'[^ATGCN]', '', sequence)
        
        return sequence
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        if not sequence:
            return 0.0
        
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def clean_sequence_input(self, sequence: str) -> str:
        """Clean sequence input - remove spaces, convert to uppercase, keep only valid bases"""
        # Remove existing spaces and convert to uppercase
        clean_sequence = ''.join(sequence.upper().split())
        
        # Only keep valid DNA bases
        valid_bases = ''.join([base for base in clean_sequence if base in 'ATGCN'])
        
        return valid_bases
    
    def format_sequence_display(self, sequence: str, line_length: int = 80) -> List[str]:
        """Format sequence for display with line numbers - no spaces, clean DNA"""
        lines = []
        
        # Calculate the maximum line number to determine padding width
        total_lines = (len(sequence) + line_length - 1) // line_length
        max_line_num = ((total_lines - 1) * line_length) + 1
        padding_width = max(5, len(str(max_line_num)))  # Minimum 5 digits, or more if needed
        
        for i in range(0, len(sequence), line_length):
            line_num = str(i + 1).zfill(padding_width)  # Zero-pad the line number
            seq_chunk = sequence[i:i + line_length]
            
            # No spaces - clean DNA sequence display
            lines.append(f"{line_num}: {seq_chunk}")
        
        return lines
    
    def format_sequence_with_reverse_complement(self, sequence: str, line_length: int = 80) -> List[str]:
        """Format sequence for display with reverse complement strand below each line"""
        lines = []
        
        # Calculate the maximum line number to determine padding width
        total_lines = (len(sequence) + line_length - 1) // line_length
        max_line_num = ((total_lines - 1) * line_length) + 1
        padding_width = max(5, len(str(max_line_num)))  # Minimum 5 digits, or more if needed
        
        # Create reverse complement
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        reverse_complement = ''.join(complement_map.get(base, base) for base in sequence[::-1])
        
        for i in range(0, len(sequence), line_length):
            line_num = str(i + 1).zfill(padding_width)  # Zero-pad the line number
            seq_chunk = sequence[i:i + line_length]
            
            # Calculate corresponding reverse complement chunk
            # For reverse complement, we need to map the position correctly
            rev_start = len(sequence) - (i + len(seq_chunk))
            rev_end = len(sequence) - i
            rev_chunk = reverse_complement[rev_start:rev_end]
            
            # Add forward strand (5' -> 3') - separate line number from sequence
            lines.append(f"{line_num}|{seq_chunk}")
            
            # Add reverse complement strand (3' -> 5') - use same line number format for alignment
            lines.append(f"{' ' * len(line_num)}|{rev_chunk[::-1]}")  # Reverse the chunk to show 3'->5' direction
            
            # Add a small separator between double strands (except for last chunk)
            if i + line_length < len(sequence):
                lines.append("")
        
        return lines
    
    def analyze_gc_content(self, sequence: str) -> Dict:
        """Detailed GC content analysis"""
        if not sequence:
            return {'overall': 0, 'distribution': []}
        
        seq_obj = Seq(sequence)
        overall_gc = gc_fraction(seq_obj) * 100
        
        # Calculate GC content in windows
        window_size = 100
        distribution = []
        
        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            if len(window) >= 10:  # Only analyze windows with sufficient length
                window_seq = Seq(window)
                window_gc = gc_fraction(window_seq) * 100
                distribution.append({
                    'position': i + 1,
                    'end_position': min(i + window_size, len(sequence)),
                    'gc_content': round(window_gc, 2)
                })
        
        return {
            'overall': round(overall_gc, 2),
            'distribution': distribution
        }
    
    def calculate_melting_temp(self, sequence: str) -> Dict:
        """Calculate melting temperature using different methods"""
        if not sequence:
            return {'basic': 0, 'salt_adjusted': 0, 'nearest_neighbor': 0}
        
        length = len(sequence)
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        
        # Basic Wallace rule: Tm = 2(A+T) + 4(G+C)
        basic_tm = 2 * at_count + 4 * gc_count
        
        # Salt-adjusted formula (more accurate for PCR conditions)
        # Tm = 81.5 + 16.6(log10[Na+]) + 0.41(%GC) - 675/length
        # Assuming 50mM Na+ and 1.5mM Mg2+ (standard PCR conditions)
        if length > 0:
            gc_percent = (gc_count / length) * 100
            salt_adjusted_tm = 81.5 + 16.6 * 1.699 + 0.41 * gc_percent - 675 / length  # log10(50mM) ≈ 1.699
        else:
            salt_adjusted_tm = 0
        
        # Nearest neighbor approximation (simplified)
        # This is a simplified version - full nearest neighbor requires dinucleotide analysis
        nn_tm = 0
        if length > 1:
            # Simplified nearest neighbor based on GC content and length
            nn_tm = 64.9 + 41 * (gc_count - 16.4) / length if length > 0 else 0
        
        return {
            'basic': round(basic_tm, 1),
            'salt_adjusted': round(salt_adjusted_tm, 1),
            'nearest_neighbor': round(nn_tm, 1) if nn_tm > 0 else round(basic_tm, 1)
        }
    
    def find_restriction_sites(self, sequence: str, enzyme_name: str = None) -> Dict:
        """Find restriction enzyme recognition sites"""
        sites = {}
        
        enzymes_to_check = {enzyme_name: self.restriction_enzymes[enzyme_name]} if enzyme_name and enzyme_name in self.restriction_enzymes else self.restriction_enzymes
        
        for enzyme, recognition_site in enzymes_to_check.items():
            positions = []
            start = 0
            
            while True:
                pos = sequence.find(recognition_site, start)
                if pos == -1:
                    break
                positions.append(pos + 1)  # 1-based indexing
                start = pos + 1
            
            # Always add enzyme data, even if no sites found
            sites[enzyme] = {
                'recognition_site': recognition_site,
                'positions': positions,
                'count': len(positions)
            }
        
        return sites
    
    def translate_reading_frames(self, sequence: str, frames: List[str]) -> Dict:
        """Translate DNA in specified reading frames"""
        translations = {}
        seq_obj = Seq(sequence)
        
        for frame in frames:
            if frame.startswith('+'):
                # Forward frames
                frame_num = int(frame[1:]) - 1
                frame_seq = seq_obj[frame_num:]
            else:
                # Reverse frames
                frame_num = int(frame[1:]) - 1
                reverse_comp = seq_obj.reverse_complement()
                frame_seq = reverse_comp[frame_num:]
            
            # Translate to protein
            protein = str(frame_seq.translate())
            
            translations[frame] = {
                'sequence': str(frame_seq),
                'protein': protein,
                'length': len(frame_seq)
            }
        
        return translations
    
    def find_amino_acids(self, sequence: str, target_aa: str) -> Dict:
        """Find specific amino acids in all reading frames"""
        results = {}
        seq_obj = Seq(sequence)
        
        # Define codon to amino acid mapping
        codon_table = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'N': ['AAT', 'AAC'], 'D': ['GAT', 'GAC'], 'C': ['TGT', 'TGC'],
            'Q': ['CAA', 'CAG'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'K': ['AAA', 'AAG'], 'M': ['ATG'], 'F': ['TTT', 'TTC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'W': ['TGG'], 'Y': ['TAT', 'TAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'STOP': ['TAA', 'TAG', 'TGA'], 'START': ['ATG']
        }
        
        if target_aa == 'ALL':
            # Find all amino acids
            for aa, codons in codon_table.items():
                if aa not in ['STOP', 'START']:
                    aa_results = self._find_aa_positions(sequence, aa, codons)
                    if aa_results['positions']:
                        results[aa] = aa_results
        else:
            # Find specific amino acid
            if target_aa in codon_table:
                codons = codon_table[target_aa]
                results[target_aa] = self._find_aa_positions(sequence, target_aa, codons)
        
        return results
    
    def _find_aa_positions(self, sequence: str, aa: str, codons: List[str]) -> Dict:
        """Helper method to find positions of specific amino acid codons"""
        positions = []
        codon_counts = {}
        
        # Convert RNA codons to DNA codons
        dna_codons = [codon.replace('U', 'T') for codon in codons]
        
        for codon in dna_codons:
            codon_positions = []
            start = 0
            while True:
                pos = sequence.find(codon, start)
                if pos == -1:
                    break
                codon_positions.append(pos + 1)  # 1-based indexing
                start = pos + 1
            
            if codon_positions:
                codon_counts[codon] = len(codon_positions)
                positions.extend([(pos, codon) for pos in codon_positions])
        
        # Sort positions
        positions.sort(key=lambda x: x[0])
        
        return {
            'amino_acid': aa,
            'positions': positions,
            'total_count': len(positions),
            'codon_usage': codon_counts
        }
    
    def get_optimal_codon_for_gc(self, amino_acid: str, target_gc: float) -> str:
        """Get optimal codon for amino acid based on GC target (copied from codon-optimizer.py)"""
        # Standard genetic code
        genetic_code = {
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
        
        # Reverse genetic code (amino acid to codons)
        aa_to_codons = {}
        for codon, aa in genetic_code.items():
            if aa not in aa_to_codons:
                aa_to_codons[aa] = []
            aa_to_codons[aa].append(codon)
        
        if amino_acid not in aa_to_codons:
            return 'NNN'
        
        available_codons = aa_to_codons[amino_acid]
        
        # Score codons by GC content preference (exact copy from codon-optimizer.py)
        scored_codons = []
        for codon in available_codons:
            codon_gc = self.calculate_gc_content(codon)
            gc_diff = abs(codon_gc - target_gc)
            # Score combines GC target preference
            score = 1 / (1 + gc_diff/10)
            scored_codons.append((codon, score))
        
        if scored_codons:
            scored_codons.sort(key=lambda x: x[1], reverse=True)
            return scored_codons[0][0]
        
        # Fallback to first available codon
        return available_codons[0]

    def optimize_gc_content(self, sequence: str, target_range: str) -> Dict:
        """Optimize GC content by rewriting the sequence (using codon-optimizer.py approach)"""
        current_gc = self.calculate_gc_content(sequence)
        
        # Parse target range and determine AGGRESSIVE target
        if target_range == 'custom':
            target_min, target_max = 45, 55  # Default balanced range
            target_gc = 50  # Balanced target
        elif '-' in target_range:
            target_min, target_max = map(int, target_range.split('-'))
            
            # AGGRESSIVE TARGETING - go for the EXTREME ends!
            if target_range == "45-55":
                # Special case: BALANCED range always targets 50% (middle)
                target_gc = 50
            elif current_gc > target_max:
                # Current GC is too high - target the LOWEST possible (minimum)
                target_gc = target_min
            elif current_gc < target_min:
                # Current GC is too low - target the HIGHEST possible (maximum)  
                target_gc = target_max
            else:
                # Already in range - target the middle (balanced)
                target_gc = (target_min + target_max) / 2
        else:
            # Handle custom slider value (direct percentage)
            try:
                target_gc = float(target_range)
                target_min = target_gc - 2.5
                target_max = target_gc + 2.5
            except:
                target_min, target_max = 45, 55
                target_gc = 50
        
        # Standard genetic code
        genetic_code = {
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
        
        # Try to optimize the sequence
        clean_seq = sequence.replace(' ', '').replace('\n', '')
        
        if len(clean_seq) % 3 != 0:
            # Not a complete coding sequence, can't optimize codons
            return {
                'current_gc': current_gc,
                'target_range': f'{target_min}-{target_max}%',
                'target_gc': target_gc,
                'optimized_sequence': None,
                'optimized_gc': None,
                'optimization_needed': current_gc < target_min or current_gc > target_max,
                'error': 'Sequence length is not divisible by 3 - cannot optimize codons'
            }
        
        # Check if optimization is actually needed (strict range check)
        if target_min <= current_gc <= target_max:
            # Already in target range, minimal optimization needed
            return {
                'current_gc': current_gc,
                'target_range': f'{target_min}-{target_max}%',
                'target_gc': target_gc,
                'optimized_sequence': clean_seq,  # Keep original
                'optimized_gc': current_gc,
                'optimization_needed': False,
                'improvement': 0,
                'changes_made': 0,
                'message': 'GC content already optimal - no changes needed'
            }
        
        # Optimize sequence using codon-optimizer.py approach
        optimized_sequence = ""
        changes_made = 0
        
        for i in range(0, len(clean_seq) - 2, 3):
            codon = clean_seq[i:i+3]
            if len(codon) == 3 and codon in genetic_code:
                aa = genetic_code[codon]
                if aa != '*':  # Don't optimize stop codons
                    optimal_codon = self.get_optimal_codon_for_gc(aa, target_gc)
                    optimized_sequence += optimal_codon
                    if optimal_codon != codon:
                        changes_made += 1
                else:
                    optimized_sequence += codon
            else:
                optimized_sequence += codon
        
        optimized_gc = self.calculate_gc_content(optimized_sequence)
        
        # Calculate improvement
        improvement = abs(target_gc - current_gc) - abs(target_gc - optimized_gc)
        
        return {
            'current_gc': current_gc,
            'target_range': f'{target_min}-{target_max}%',
            'target_gc': target_gc,
            'optimized_sequence': optimized_sequence,
            'optimized_gc': optimized_gc,
            'optimization_needed': current_gc < target_min or current_gc > target_max,
            'improvement': improvement,
            'changes_made': changes_made
        }
    
    def optimize_codon_usage(self, sequence: str, organism: str) -> Dict:
        """Optimize codon usage for specific organism - ACTUALLY REPLACES CODONS!"""
        # Import the comprehensive codon usage tables from codon-optimizer.py
        import sys
        import os
        import importlib.util
        
        # Load the codon-optimizer.py module
        codon_optimizer_path = os.path.join(os.path.dirname(__file__), 'codon-optimizer.py')
        spec = importlib.util.spec_from_file_location("codon_optimizer", codon_optimizer_path)
        codon_optimizer_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(codon_optimizer_module)
        CodonOptimizer = codon_optimizer_module.CodonOptimizer
        
        try:
            # Create codon optimizer instance
            optimizer = CodonOptimizer()
            
            # Check if organism is supported
            if organism not in optimizer.codon_usage:
                available_organisms = list(optimizer.codon_usage.keys())
                return {
                    'organism': organism,
                    'error': f'Organism "{organism}" not supported. Available: {", ".join(available_organisms)}',
                    'success': False
                }
            
            # Clean sequence
            clean_seq = sequence.replace(' ', '').replace('\n', '').upper()
            
            # Calculate original protein sequence to verify it stays the same
            original_protein = optimizer.translate_dna(clean_seq)
            original_gc = optimizer.calculate_gc_content(clean_seq)
            
            # Optimize the sequence using the professional codon optimizer
            optimized_sequence = ""
            changes_made = 0
            codon_changes = []
            
            for i in range(0, len(clean_seq) - 2, 3):
                codon = clean_seq[i:i+3]
                if len(codon) == 3 and codon in optimizer.genetic_code:
                    aa = optimizer.genetic_code[codon]
                    if aa != '*':  # Don't optimize stop codons
                        # Get optimal codon for this amino acid in the target organism
                        optimal_codon = optimizer.get_optimal_codon(aa, organism)
                        optimized_sequence += optimal_codon
                        
                        # Track changes
                        if optimal_codon != codon:
                            changes_made += 1
                            codon_changes.append({
                                'position': i + 1,
                                'amino_acid': aa,
                                'original_codon': codon,
                                'optimized_codon': optimal_codon,
                                'improvement': 'Organism-preferred codon'
                            })
                    else:
                        optimized_sequence += codon  # Keep stop codons as-is
                else:
                    optimized_sequence += codon  # Keep incomplete codons as-is
            
            # Verify protein sequence is identical
            optimized_protein = optimizer.translate_dna(optimized_sequence)
            optimized_gc = optimizer.calculate_gc_content(optimized_sequence)
            
            # Calculate optimization metrics
            total_codons = len([codon for codon in [clean_seq[i:i+3] for i in range(0, len(clean_seq)-2, 3)] 
                              if len(codon) == 3 and codon in optimizer.genetic_code and optimizer.genetic_code[codon] != '*'])
            
            optimization_percentage = (changes_made / total_codons * 100) if total_codons > 0 else 0
            
            return {
                'success': True,
                'organism': organism,
                'original_sequence': clean_seq,
                'optimized_sequence': optimized_sequence,
                'original_protein': original_protein,
                'optimized_protein': optimized_protein,
                'protein_identical': original_protein == optimized_protein,
                'original_gc': round(original_gc, 1),
                'optimized_gc': round(optimized_gc, 1),
                'gc_change': round(optimized_gc - original_gc, 1),
                'total_codons': total_codons,
                'changes_made': changes_made,
                'optimization_percentage': round(optimization_percentage, 1),
                'codon_changes': codon_changes[:10],  # Show first 10 changes
                'message': f'Optimized {changes_made} codons ({optimization_percentage:.1f}%) for {organism} expression'
            }
            
        except ImportError:
            # Fallback to simplified codon tables if codon-optimizer.py not available
            return self._optimize_codon_usage_fallback(sequence, organism)
        except Exception as e:
            return {
                'success': False,
                'organism': organism,
                'error': f'Codon optimization failed: {str(e)}'
            }
    
    def _optimize_codon_usage_fallback(self, sequence: str, organism: str) -> Dict:
        """Fallback codon optimization using simplified tables"""
        return {
            'success': False,
            'organism': organism,
            'error': 'Codon optimization requires the full codon-optimizer.py module. Please ensure it is available.',
            'message': 'Codon optimization is currently unavailable'
        }
    
    def _generate_codon_recommendations(self, current_codons: List[str], optimal_codons: List[str], organism: str) -> List[Dict]:
        """Generate codon optimization recommendations"""
        recommendations = []
        
        # Count suboptimal codons
        suboptimal_count = 0
        codon_changes = {}
        
        for i, (current, optimal) in enumerate(zip(current_codons, optimal_codons)):
            if current != optimal:
                suboptimal_count += 1
                change_key = f"{current}→{optimal}"
                if change_key not in codon_changes:
                    codon_changes[change_key] = 0
                codon_changes[change_key] += 1
        
        if suboptimal_count == 0:
            recommendations.append({
                'type': 'success',
                'message': f'Sequence is already optimized for {organism}',
                'priority': 'low'
            })
        else:
            recommendations.append({
                'type': 'optimization',
                'message': f'{suboptimal_count} codons could be optimized for {organism}',
                'priority': 'high' if suboptimal_count > len(current_codons) * 0.3 else 'medium',
                'suggested_changes': dict(list(codon_changes.items())[:5])  # Top 5 changes
            })
        
        return recommendations
    
    def apply_color_scheme(self, formatted_lines: List[str], color_scheme: str) -> List[str]:
        """Apply color scheme to formatted sequence lines"""
        if color_scheme not in self.color_schemes:
            return formatted_lines
        
        colors = self.color_schemes[color_scheme]
        colored_lines = []
        
        for line in formatted_lines:
            # Split line number from sequence
            parts = line.split(': ', 1)
            if len(parts) == 2:
                line_num, sequence_part = parts
                
                # Apply colors to nucleotides
                colored_sequence = sequence_part
                for nucleotide, color_class in colors.items():
                    colored_sequence = colored_sequence.replace(
                        nucleotide, 
                        f'<span class="{color_class}">{nucleotide}</span>'
                    )
                
                colored_lines.append(f"{line_num}: {colored_sequence}")
            else:
                colored_lines.append(line)
        
        return colored_lines
    
    def apply_amino_acid_highlighting(self, formatted_lines: List[str], sequence: str, target_amino_acids: str) -> List[str]:
        """Apply highlighting to codons that code for target amino acids"""
        # Standard genetic code
        genetic_code = {
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
        
        # Get target amino acids
        if target_amino_acids == 'ALL':
            target_aas = set(genetic_code.values()) - {'*'}  # All except stop codons
        elif target_amino_acids == 'START':
            target_aas = {'M'}  # Start codon
        elif target_amino_acids == 'STOP':
            target_aas = {'*'}  # Stop codons
        else:
            target_aas = {target_amino_acids}
        
        # Find codons that code for target amino acids
        target_codons = []
        for codon, aa in genetic_code.items():
            if aa in target_aas:
                target_codons.append(codon)
        
        # Find positions of target codons in the sequence
        highlight_positions = set()
        clean_sequence = sequence.replace(' ', '').replace('\n', '')
        
        for i in range(0, len(clean_sequence) - 2, 3):
            codon = clean_sequence[i:i+3]
            if len(codon) == 3 and codon in target_codons:
                # Mark all three positions of this codon for highlighting
                highlight_positions.update([i, i+1, i+2])
        
        # Apply highlighting to formatted lines
        highlighted_lines = []
        sequence_position = 0
        
        for line in formatted_lines:
            parts = line.split(': ', 1)
            if len(parts) == 2:
                line_num, sequence_part = parts
                highlighted_sequence = ""
                
                # Process each character in the sequence part
                for char in sequence_part:
                    if char in 'ATGCN':
                        if sequence_position in highlight_positions:
                            # Highlight this nucleotide
                            highlighted_sequence += f'<span class="bg-yellow-300 text-black font-bold px-0.5 rounded">{char}</span>'
                        else:
                            highlighted_sequence += char
                        sequence_position += 1
                    else:
                        # Space or other character
                        highlighted_sequence += char
                
                highlighted_lines.append(f"{line_num}: {highlighted_sequence}")
            else:
                highlighted_lines.append(line)
        
        return highlighted_lines
    
    def handle_special_restriction_searches(self, sequence: str, search_type: str) -> Dict:
        """Handle special restriction enzyme search categories"""
        results = {}
        
        if search_type == 'ALL_COMMON':
            # Common 6-base cutters
            common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XbaI', 'SalI', 'PstI', 'SmaI', 'KpnI']
            for enzyme in common_enzymes:
                if enzyme in self.restriction_enzymes:
                    sites = self._find_single_enzyme_sites(sequence, enzyme, self.restriction_enzymes[enzyme])
                    results[enzyme] = sites
        
        elif search_type == 'ALL_RARE':
            # 8-base rare cutters
            rare_enzymes = ['NotI', 'AsiSI', 'FseI', 'PacI', 'SbfI', 'SwaI', 'PmeI']
            for enzyme in rare_enzymes:
                if enzyme in self.restriction_enzymes:
                    sites = self._find_single_enzyme_sites(sequence, enzyme, self.restriction_enzymes[enzyme])
                    results[enzyme] = sites
        
        elif search_type == 'ALL_BLUNT':
            # Blunt end cutters
            blunt_enzymes = ['EcoRV', 'SmaI', 'HpaI', 'ScaI', 'StuI', 'NruI', 'SspI']
            for enzyme in blunt_enzymes:
                if enzyme in self.restriction_enzymes:
                    sites = self._find_single_enzyme_sites(sequence, enzyme, self.restriction_enzymes[enzyme])
                    results[enzyme] = sites
        
        elif search_type == 'ALL_STICKY':
            # Sticky end cutters (most 6-base cutters)
            sticky_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XbaI', 'SalI', 'PstI', 'KpnI', 'XhoI']
            for enzyme in sticky_enzymes:
                if enzyme in self.restriction_enzymes:
                    sites = self._find_single_enzyme_sites(sequence, enzyme, self.restriction_enzymes[enzyme])
                    results[enzyme] = sites
        
        return results
    
    def _find_single_enzyme_sites(self, sequence: str, enzyme_name: str, recognition_site: str) -> Dict:
        """Find sites for a single restriction enzyme"""
        positions = []
        start = 0
        
        while True:
            pos = sequence.find(recognition_site, start)
            if pos == -1:
                break
            positions.append(pos + 1)  # 1-based indexing
            start = pos + 1
        
        return {
            'recognition_site': recognition_site,
            'positions': positions,
            'count': len(positions)
        }
    
    def analyze_sequence(self, sequence: str, options: Dict) -> Dict:
        """Main sequence analysis function"""
        try:
            # Clean the sequence
            clean_seq = self.clean_sequence(sequence)
            
            if not clean_seq:
                return {
                    'success': False,
                    'error': 'Invalid or empty DNA sequence'
                }
            
            # STEP 1: Codon Usage Optimization FIRST (if requested)
            working_sequence = clean_seq
            codon_optimization_result = None
            
            if options.get('codon_optimizer'):
                organism = options['codon_optimizer']
                codon_optimization_result = self.optimize_codon_usage(clean_seq, organism)
                
                # Use codon-optimized sequence if available
                if codon_optimization_result.get('success') and codon_optimization_result.get('optimized_sequence'):
                    working_sequence = codon_optimization_result['optimized_sequence']
            
            # STEP 2: GC Content Optimization SECOND (using codon-optimized sequence)
            gc_optimization_result = None
            
            if options.get('gc_optimizer'):
                target_range = options['gc_optimizer']
                # Apply GC optimization to the codon-optimized sequence (or original if no codon opt)
                gc_optimization_result = self.optimize_gc_content(working_sequence, target_range)
                
                # Use GC-optimized sequence if available
                if gc_optimization_result.get('optimized_sequence'):
                    working_sequence = gc_optimization_result['optimized_sequence']
            
            # Basic sequence info (using working sequence)
            sequence_info = {
                'length': len(working_sequence),
                'composition': {
                    'A': working_sequence.count('A'),
                    'T': working_sequence.count('T'),
                    'G': working_sequence.count('G'),
                    'C': working_sequence.count('C'),
                    'N': working_sequence.count('N')
                }
            }
            
            # Format sequence for display (using working sequence) - EXTENDED WIDTH, Clean DNA
            if options.get('show_reverse_complement', False):
                formatted_lines = self.format_sequence_with_reverse_complement(working_sequence, line_length=190)
            else:
                formatted_lines = self.format_sequence_display(working_sequence, line_length=190)
            
            # GC content analysis (using working sequence)
            gc_analysis = self.analyze_gc_content(working_sequence)
            sequence_info['gc_content'] = gc_analysis['overall']
            
            # Melting temperature calculation (always calculated for display)
            melting_temp = self.calculate_melting_temp(working_sequence)
            sequence_info['melting_temp'] = melting_temp
            
            results = {
                'sequence_info': sequence_info,
                'formatted_sequence': formatted_lines,
                'gc_analysis': gc_analysis,
                'melting_temp': melting_temp
            }
            
            # Optional analyses based on options
            if options.get('analyze_gc_content', False):
                results['detailed_gc'] = gc_analysis
            
            # Restriction enzyme analysis (using working sequence)
            if options.get('restriction_enzyme'):
                enzyme = options['restriction_enzyme']
                if enzyme.startswith('ALL_'):
                    # Handle special search categories
                    results['restriction_sites'] = self.handle_special_restriction_searches(working_sequence, enzyme)
                else:
                    # Handle single enzyme
                    results['restriction_sites'] = self.find_restriction_sites(working_sequence, enzyme)
            
            # Reading frame translations (using working sequence)
            if options.get('reading_frames'):
                frames = options['reading_frames']
                results['translations'] = self.translate_reading_frames(working_sequence, frames)
            
            # Process sequence display with highlighting and coloring
            display_sequence = formatted_lines
            
            # Apply amino acid highlighting first (if selected)
            if options.get('amino_acids_finder'):
                target_aa = options['amino_acids_finder']
                results['amino_acid_analysis'] = self.find_amino_acids(working_sequence, target_aa)
                # Apply highlighting to the sequence
                display_sequence = self.apply_amino_acid_highlighting(display_sequence, working_sequence, target_aa)
                results['amino_acid_highlighted'] = target_aa
            
            # Apply color scheme on top of highlighting (if selected)
            if options.get('color_scheme'):
                color_scheme = options['color_scheme']
                display_sequence = self.apply_color_scheme(display_sequence, color_scheme)
                results['color_scheme_used'] = color_scheme
            
            # Update the sequence display
            if display_sequence != formatted_lines:
                results['colored_sequence'] = display_sequence
            
            # Add optimization results if they were performed
            if codon_optimization_result:
                results['codon_optimization'] = codon_optimization_result
                
            if gc_optimization_result:
                results['gc_optimization'] = gc_optimization_result
                
            # Add info about sequence optimization
            if codon_optimization_result or gc_optimization_result:
                results['sequence_optimized'] = True
                results['original_sequence_info'] = {
                    'length': len(clean_seq),
                    'gc_content': self.calculate_gc_content(clean_seq)
                }
                
                # Show optimization pipeline
                optimization_steps = []
                if codon_optimization_result and codon_optimization_result.get('success'):
                    optimization_steps.append(f"Codon optimized for {codon_optimization_result.get('organism', 'unknown')}")
                if gc_optimization_result and gc_optimization_result.get('optimized_sequence'):
                    optimization_steps.append(f"GC optimized to {gc_optimization_result.get('target_range', 'unknown')}")
                
                results['optimization_pipeline'] = optimization_steps
            
            return {
                'success': True,
                'original_sequence': sequence,
                'clean_sequence': clean_seq,
                'results': results
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }


def main():
    """Main function for command-line and stdin execution"""
    sequencer = DNASequencer()
    
    try:
        # Check if we have command line arguments (legacy mode)
        if len(sys.argv) > 1:
            # Legacy command-line mode
            parser = argparse.ArgumentParser(description='DNA Sequencer - Advanced sequence analysis')
            parser.add_argument('--sequence', required=True, help='DNA sequence to analyze')
            parser.add_argument('--file-content', help='File content for file-based analysis')
            parser.add_argument('--file-format', default='txt', help='File format (fasta, genbank, txt)')
            parser.add_argument('--options', default='{}', help='Analysis options as JSON string')
            
            args = parser.parse_args()
            
            # Parse options
            options = json.loads(args.options) if args.options else {}
            
            # Handle file content if provided
            if args.file_content:
                sequences = sequencer.parse_file_content(args.file_content, args.file_format)
                if not sequences:
                    raise ValueError("No valid sequences found in file")
                
                # Use first sequence for analysis
                sequence = sequences[0]['sequence']
                sequence_id = sequences[0]['id']
                sequence_description = sequences[0]['description']
                file_format = args.file_format
            else:
                sequence = args.sequence
                # Clean manual sequence input - remove spaces, keep only valid bases
                sequence = sequencer.clean_sequence_input(sequence)
                sequence_id = 'manual_input'
                sequence_description = 'Manually entered sequence'
                file_format = None
        else:
            # New stdin mode - read JSON input from stdin
            stdin_data = sys.stdin.read().strip()
            if not stdin_data:
                raise ValueError("No input data provided")
            
            input_data = json.loads(stdin_data)
            
            # Extract data from JSON input
            sequence = input_data.get('sequence', '')
            file_content = input_data.get('fileContent')
            file_format = input_data.get('fileFormat', 'txt')
            options = input_data.get('options', {})
            
            # Clean manual sequence input - remove spaces, keep only valid bases
            if sequence and not file_content:
                sequence = sequencer.clean_sequence_input(sequence)
            
            # Handle file content if provided
            if file_content:
                sequences = sequencer.parse_file_content(file_content, file_format)
                if not sequences:
                    raise ValueError("No valid sequences found in file")
                
                # Use first sequence for analysis
                sequence = sequences[0]['sequence']
                sequence_id = sequences[0]['id']
                sequence_description = sequences[0]['description']
            else:
                sequence_id = 'manual_input'
                sequence_description = 'Manually entered sequence'
                file_format = None
        
        # Perform analysis
        result = sequencer.analyze_sequence(sequence, options)
        
        if result['success']:
            # Add file info if from file
            if file_format:
                result['file_info'] = {
                    'id': sequence_id,
                    'description': sequence_description,
                    'format': file_format
                }
        
        print(json.dumps(result))
        
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