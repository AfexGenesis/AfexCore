#!/usr/bin/env python3
"""
Codon Optimizer - Advanced DNA sequence optimization tool
Optimizes codon usage for different organisms and GC content targets
"""

import argparse
import json
import sys
import re
from typing import Dict, List, Tuple, Optional
from collections import Counter
import random

class CodonOptimizer:
    def __init__(self):
        """Initialize the Codon Optimizer with codon tables and organism preferences."""
        
        # Standard genetic code
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
        
        # Reverse genetic code (amino acid to codons)
        self.aa_to_codons = {}
        for codon, aa in self.genetic_code.items():
            if aa not in self.aa_to_codons:
                self.aa_to_codons[aa] = []
            self.aa_to_codons[aa].append(codon)
        
        # Organism-specific codon usage preferences (relative frequencies)
        self.codon_usage = {
            'e-coli': {
                'F': {'TTT': 0.58, 'TTC': 0.42},
                'L': {'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47},
                'S': {'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25},
                'Y': {'TAT': 0.59, 'TAC': 0.41},
                'C': {'TGT': 0.46, 'TGC': 0.54},
                'W': {'TGG': 1.00},
                'P': {'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49},
                'H': {'CAT': 0.57, 'CAC': 0.43},
                'Q': {'CAA': 0.34, 'CAG': 0.66},
                'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04},
                'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11},
                'M': {'ATG': 1.00},
                'T': {'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25},
                'N': {'AAT': 0.49, 'AAC': 0.51},
                'K': {'AAA': 0.74, 'AAG': 0.26},
                'V': {'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35},
                'A': {'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33},
                'D': {'GAT': 0.63, 'GAC': 0.37},
                'E': {'GAA': 0.68, 'GAG': 0.32},
                'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15},
                '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30}
            },
            'human': {
                'F': {'TTT': 0.45, 'TTC': 0.55},
                'L': {'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41},
                'S': {'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.05, 'AGT': 0.15, 'AGC': 0.24},
                'Y': {'TAT': 0.43, 'TAC': 0.57},
                'C': {'TGT': 0.45, 'TGC': 0.55},
                'W': {'TGG': 1.00},
                'P': {'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11},
                'H': {'CAT': 0.41, 'CAC': 0.59},
                'Q': {'CAA': 0.25, 'CAG': 0.75},
                'R': {'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.20},
                'I': {'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16},
                'M': {'ATG': 1.00},
                'T': {'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12},
                'N': {'AAT': 0.46, 'AAC': 0.54},
                'K': {'AAA': 0.42, 'AAG': 0.58},
                'V': {'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47},
                'A': {'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11},
                'D': {'GAT': 0.46, 'GAC': 0.54},
                'E': {'GAA': 0.42, 'GAG': 0.58},
                'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25},
                '*': {'TAA': 0.28, 'TAG': 0.20, 'TGA': 0.52}
            },
            'yeast': {
                'F': {'TTT': 0.59, 'TTC': 0.41},
                'L': {'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06, 'CTA': 0.14, 'CTG': 0.11},
                'S': {'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.10, 'AGT': 0.16, 'AGC': 0.11},
                'Y': {'TAT': 0.56, 'TAC': 0.44},
                'C': {'TGT': 0.63, 'TGC': 0.37},
                'W': {'TGG': 1.00},
                'P': {'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.42, 'CCG': 0.12},
                'H': {'CAT': 0.64, 'CAC': 0.36},
                'Q': {'CAA': 0.69, 'CAG': 0.31},
                'R': {'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21},
                'I': {'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27},
                'M': {'ATG': 1.00},
                'T': {'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.30, 'ACG': 0.14},
                'N': {'AAT': 0.59, 'AAC': 0.41},
                'K': {'AAA': 0.58, 'AAG': 0.42},
                'V': {'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21, 'GTG': 0.19},
                'A': {'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11},
                'D': {'GAT': 0.65, 'GAC': 0.35},
                'E': {'GAA': 0.70, 'GAG': 0.30},
                'G': {'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12},
                '*': {'TAA': 0.47, 'TAG': 0.23, 'TGA': 0.30}
            }
        }
        
        # Add more organisms with simplified preferences
        self.codon_usage.update({
            'arabidopsis': self._create_plant_codon_usage(),
            'mouse': self.codon_usage['human'],  # Similar to human
            'cho-cells': self.codon_usage['human'],  # CHO cells similar to human
            'zebrafish': self._create_vertebrate_codon_usage(),
            'c-elegans': self._create_nematode_codon_usage(),
            'baculovirus': self._create_insect_codon_usage()
        })
        
        # Common restriction enzyme sites to avoid
        self.restriction_sites = [
            'GAATTC',  # EcoRI
            'GGATCC',  # BamHI
            'AAGCTT',  # HindIII
            'CTGCAG',  # PstI
            'GTCGAC',  # SalI
            'CCCGGG',  # SmaI
            'TCTAGA',  # XbaI
            'TCTAGG',  # XbaI (partial)
            'GCGGCCGC',  # NotI
            'GGCCGGCC'   # NotI (partial)
        ]
    
    def _create_plant_codon_usage(self) -> Dict:
        """Create plant-specific codon usage (Arabidopsis-like)."""
        return {
            'F': {'TTT': 0.52, 'TTC': 0.48},
            'L': {'TTA': 0.08, 'TTG': 0.15, 'CTT': 0.18, 'CTC': 0.12, 'CTA': 0.08, 'CTG': 0.39},
            'S': {'TCT': 0.20, 'TCC': 0.18, 'TCA': 0.16, 'TCG': 0.12, 'AGT': 0.18, 'AGC': 0.16},
            'Y': {'TAT': 0.48, 'TAC': 0.52},
            'C': {'TGT': 0.48, 'TGC': 0.52},
            'W': {'TGG': 1.00},
            'P': {'CCT': 0.25, 'CCC': 0.20, 'CCA': 0.30, 'CCG': 0.25},
            'H': {'CAT': 0.48, 'CAC': 0.52},
            'Q': {'CAA': 0.35, 'CAG': 0.65},
            'R': {'CGT': 0.12, 'CGC': 0.15, 'CGA': 0.08, 'CGG': 0.15, 'AGA': 0.25, 'AGG': 0.25},
            'I': {'ATT': 0.40, 'ATC': 0.35, 'ATA': 0.25},
            'M': {'ATG': 1.00},
            'T': {'ACT': 0.25, 'ACC': 0.30, 'ACA': 0.25, 'ACG': 0.20},
            'N': {'AAT': 0.48, 'AAC': 0.52},
            'K': {'AAA': 0.45, 'AAG': 0.55},
            'V': {'GTT': 0.22, 'GTC': 0.18, 'GTA': 0.15, 'GTG': 0.45},
            'A': {'GCT': 0.28, 'GCC': 0.32, 'GCA': 0.25, 'GCG': 0.15},
            'D': {'GAT': 0.52, 'GAC': 0.48},
            'E': {'GAA': 0.48, 'GAG': 0.52},
            'G': {'GGT': 0.25, 'GGC': 0.30, 'GGA': 0.25, 'GGG': 0.20},
            '*': {'TAA': 0.40, 'TAG': 0.25, 'TGA': 0.35}
        }
    
    def _create_vertebrate_codon_usage(self) -> Dict:
        """Create vertebrate-specific codon usage (zebrafish-like)."""
        return {
            'F': {'TTT': 0.42, 'TTC': 0.58},
            'L': {'TTA': 0.06, 'TTG': 0.12, 'CTT': 0.12, 'CTC': 0.22, 'CTA': 0.06, 'CTG': 0.42},
            'S': {'TCT': 0.16, 'TCC': 0.24, 'TCA': 0.14, 'TCG': 0.06, 'AGT': 0.14, 'AGC': 0.26},
            'Y': {'TAT': 0.40, 'TAC': 0.60},
            'C': {'TGT': 0.42, 'TGC': 0.58},
            'W': {'TGG': 1.00},
            'P': {'CCT': 0.26, 'CCC': 0.36, 'CCA': 0.26, 'CCG': 0.12},
            'H': {'CAT': 0.38, 'CAC': 0.62},
            'Q': {'CAA': 0.22, 'CAG': 0.78},
            'R': {'CGT': 0.06, 'CGC': 0.20, 'CGA': 0.10, 'CGG': 0.22, 'AGA': 0.18, 'AGG': 0.24},
            'I': {'ATT': 0.34, 'ATC': 0.50, 'ATA': 0.16},
            'M': {'ATG': 1.00},
            'T': {'ACT': 0.22, 'ACC': 0.38, 'ACA': 0.26, 'ACG': 0.14},
            'N': {'AAT': 0.44, 'AAC': 0.56},
            'K': {'AAA': 0.38, 'AAG': 0.62},
            'V': {'GTT': 0.16, 'GTC': 0.26, 'GTA': 0.10, 'GTG': 0.48},
            'A': {'GCT': 0.24, 'GCC': 0.42, 'GCA': 0.22, 'GCG': 0.12},
            'D': {'GAT': 0.44, 'GAC': 0.56},
            'E': {'GAA': 0.40, 'GAG': 0.60},
            'G': {'GGT': 0.14, 'GGC': 0.36, 'GGA': 0.24, 'GGG': 0.26},
            '*': {'TAA': 0.26, 'TAG': 0.18, 'TGA': 0.56}
        }
    
    def _create_nematode_codon_usage(self) -> Dict:
        """Create nematode-specific codon usage (C. elegans-like)."""
        return {
            'F': {'TTT': 0.58, 'TTC': 0.42},
            'L': {'TTA': 0.12, 'TTG': 0.18, 'CTT': 0.15, 'CTC': 0.10, 'CTA': 0.08, 'CTG': 0.37},
            'S': {'TCT': 0.22, 'TCC': 0.14, 'TCA': 0.18, 'TCG': 0.08, 'AGT': 0.20, 'AGC': 0.18},
            'Y': {'TAT': 0.55, 'TAC': 0.45},
            'C': {'TGT': 0.55, 'TGC': 0.45},
            'W': {'TGG': 1.00},
            'P': {'CCT': 0.28, 'CCC': 0.18, 'CCA': 0.35, 'CCG': 0.19},
            'H': {'CAT': 0.55, 'CAC': 0.45},
            'Q': {'CAA': 0.45, 'CAG': 0.55},
            'R': {'CGT': 0.18, 'CGC': 0.12, 'CGA': 0.12, 'CGG': 0.08, 'AGA': 0.30, 'AGG': 0.20},
            'I': {'ATT': 0.45, 'ATC': 0.32, 'ATA': 0.23},
            'M': {'ATG': 1.00},
            'T': {'ACT': 0.28, 'ACC': 0.25, 'ACA': 0.30, 'ACG': 0.17},
            'N': {'AAT': 0.52, 'AAC': 0.48},
            'K': {'AAA': 0.55, 'AAG': 0.45},
            'V': {'GTT': 0.32, 'GTC': 0.18, 'GTA': 0.18, 'GTG': 0.32},
            'A': {'GCT': 0.32, 'GCC': 0.25, 'GCA': 0.28, 'GCG': 0.15},
            'D': {'GAT': 0.58, 'GAC': 0.42},
            'E': {'GAA': 0.58, 'GAG': 0.42},
            'G': {'GGT': 0.32, 'GGC': 0.22, 'GGA': 0.28, 'GGG': 0.18},
            '*': {'TAA': 0.45, 'TAG': 0.25, 'TGA': 0.30}
        }
    
    def _create_insect_codon_usage(self) -> Dict:
        """Create insect cell codon usage (baculovirus system)."""
        return {
            'F': {'TTT': 0.48, 'TTC': 0.52},
            'L': {'TTA': 0.10, 'TTG': 0.14, 'CTT': 0.16, 'CTC': 0.18, 'CTA': 0.08, 'CTG': 0.34},
            'S': {'TCT': 0.18, 'TCC': 0.20, 'TCA': 0.16, 'TCG': 0.08, 'AGT': 0.18, 'AGC': 0.20},
            'Y': {'TAT': 0.46, 'TAC': 0.54},
            'C': {'TGT': 0.46, 'TGC': 0.54},
            'W': {'TGG': 1.00},
            'P': {'CCT': 0.24, 'CCC': 0.26, 'CCA': 0.28, 'CCG': 0.22},
            'H': {'CAT': 0.44, 'CAC': 0.56},
            'Q': {'CAA': 0.32, 'CAG': 0.68},
            'R': {'CGT': 0.10, 'CGC': 0.16, 'CGA': 0.12, 'CGG': 0.18, 'AGA': 0.22, 'AGG': 0.22},
            'I': {'ATT': 0.38, 'ATC': 0.42, 'ATA': 0.20},
            'M': {'ATG': 1.00},
            'T': {'ACT': 0.24, 'ACC': 0.32, 'ACA': 0.26, 'ACG': 0.18},
            'N': {'AAT': 0.46, 'AAC': 0.54},
            'K': {'AAA': 0.44, 'AAG': 0.56},
            'V': {'GTT': 0.20, 'GTC': 0.22, 'GTA': 0.14, 'GTG': 0.44},
            'A': {'GCT': 0.26, 'GCC': 0.34, 'GCA': 0.24, 'GCG': 0.16},
            'D': {'GAT': 0.48, 'GAC': 0.52},
            'E': {'GAA': 0.46, 'GAG': 0.54},
            'G': {'GGT': 0.20, 'GGC': 0.32, 'GGA': 0.26, 'GGG': 0.22},
            '*': {'TAA': 0.35, 'TAG': 0.22, 'TGA': 0.43}
        }
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and validate DNA sequence."""
        # Remove whitespace and convert to uppercase
        sequence = re.sub(r'\s+', '', sequence.upper())
        
        # Remove non-DNA characters
        sequence = re.sub(r'[^ATGC]', '', sequence)
        
        return sequence
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        if not sequence:
            return 0.0
        
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def translate_dna(self, sequence: str) -> str:
        """Translate DNA sequence to amino acids."""
        protein = ""
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                aa = self.genetic_code.get(codon, 'X')
                if aa == '*':
                    break
                protein += aa
        return protein
    
    def get_optimal_codon(self, amino_acid: str, organism: str, gc_target: Optional[float] = None) -> str:
        """Get optimal codon for amino acid based on organism and GC target."""
        if amino_acid not in self.aa_to_codons:
            return 'NNN'
        
        available_codons = self.aa_to_codons[amino_acid]
        
        if organism in self.codon_usage and amino_acid in self.codon_usage[organism]:
            # Use organism-specific preferences
            codon_prefs = self.codon_usage[organism][amino_acid]
            
            if gc_target is not None:
                # Filter codons by GC content preference
                scored_codons = []
                for codon in available_codons:
                    if codon in codon_prefs:
                        codon_gc = self.calculate_gc_content(codon)
                        gc_diff = abs(codon_gc - gc_target)
                        # Score combines organism preference and GC target
                        score = codon_prefs[codon] * (1 / (1 + gc_diff/10))
                        scored_codons.append((codon, score))
                
                if scored_codons:
                    scored_codons.sort(key=lambda x: x[1], reverse=True)
                    return scored_codons[0][0]
            
            # Use organism preference only
            best_codon = max(codon_prefs.keys(), key=lambda x: codon_prefs[x])
            return best_codon
        
        # Fallback to first available codon
        return available_codons[0]
    
    def has_restriction_sites(self, sequence: str) -> List[str]:
        """Check for restriction enzyme sites in sequence."""
        found_sites = []
        for site in self.restriction_sites:
            if site in sequence:
                found_sites.append(site)
        return found_sites
    
    def remove_restriction_sites(self, sequence: str, organism: str) -> str:
        """Remove restriction sites by changing codons."""
        modified_sequence = sequence
        changes_made = 0
        
        for site in self.restriction_sites:
            while site in modified_sequence:
                # Find the position of the restriction site
                pos = modified_sequence.find(site)
                
                # Try to modify codons overlapping with the site
                for codon_start in range(max(0, pos - 2), min(len(modified_sequence) - 2, pos + len(site)), 3):
                    if codon_start % 3 == 0:  # Ensure we're at a codon boundary
                        codon = modified_sequence[codon_start:codon_start + 3]
                        if len(codon) == 3 and codon in self.genetic_code:
                            aa = self.genetic_code[codon]
                            if aa != '*':  # Don't change stop codons
                                # Get alternative codon
                                alternatives = [c for c in self.aa_to_codons[aa] if c != codon]
                                if alternatives:
                                    new_codon = self.get_optimal_codon(aa, organism)
                                    if new_codon != codon:
                                        modified_sequence = (modified_sequence[:codon_start] + 
                                                           new_codon + 
                                                           modified_sequence[codon_start + 3:])
                                        changes_made += 1
                                        break
                
                # If we couldn't remove the site, break to avoid infinite loop
                if site in modified_sequence and changes_made == 0:
                    break
        
        return modified_sequence
    
    def optimize_gc_content(self, sequence: str, target_gc: float, organism: str) -> Tuple[str, int]:
        """Optimize sequence for target GC content."""
        optimized_sequence = ""
        changes_made = 0
        
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3 and codon in self.genetic_code:
                aa = self.genetic_code[codon]
                if aa != '*':  # Don't optimize stop codons
                    optimal_codon = self.get_optimal_codon(aa, organism, target_gc)
                    optimized_sequence += optimal_codon
                    if optimal_codon != codon:
                        changes_made += 1
                else:
                    optimized_sequence += codon
            else:
                optimized_sequence += codon
        
        return optimized_sequence, changes_made
    
    def optimize_codon_usage(self, sequence: str, organism: str) -> Tuple[str, int]:
        """Optimize sequence for organism-specific codon usage."""
        optimized_sequence = ""
        changes_made = 0
        
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3 and codon in self.genetic_code:
                aa = self.genetic_code[codon]
                if aa != '*':  # Don't optimize stop codons
                    optimal_codon = self.get_optimal_codon(aa, organism)
                    optimized_sequence += optimal_codon
                    if optimal_codon != codon:
                        changes_made += 1
                else:
                    optimized_sequence += codon
            else:
                optimized_sequence += codon
        
        return optimized_sequence, changes_made
    
    def calculate_cai(self, sequence: str, organism: str) -> float:
        """Calculate Codon Adaptation Index (CAI)."""
        if organism not in self.codon_usage:
            return 0.0
        
        cai_values = []
        
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3 and codon in self.genetic_code:
                aa = self.genetic_code[codon]
                if aa != '*' and aa in self.codon_usage[organism]:
                    codon_freq = self.codon_usage[organism][aa].get(codon, 0.1)
                    max_freq = max(self.codon_usage[organism][aa].values())
                    relative_adaptiveness = codon_freq / max_freq if max_freq > 0 else 0.1
                    cai_values.append(relative_adaptiveness)
        
        if not cai_values:
            return 0.0
        
        # Geometric mean
        product = 1.0
        for val in cai_values:
            product *= max(val, 0.01)  # Avoid zero values
        
        return product ** (1.0 / len(cai_values))
    
    def optimize_sequence(self, sequence: str, gc_target: Optional[str] = None, 
                         organism: Optional[str] = None, avoid_rare_codons: bool = True,
                         remove_restriction: bool = False, optimize_folding: bool = False) -> Dict:
        """Main optimization function."""
        
        # Clean input sequence
        clean_seq = self.clean_sequence(sequence)
        
        if not clean_seq:
            return {
                'success': False,
                'error': 'Invalid or empty DNA sequence'
            }
        
        if len(clean_seq) % 3 != 0:
            return {
                'success': False,
                'error': 'DNA sequence length must be divisible by 3 (complete codons)'
            }
        
        # Translate to check if it's a valid coding sequence
        protein = self.translate_dna(clean_seq)
        if not protein:
            return {
                'success': False,
                'error': 'Sequence does not translate to a valid protein'
            }
        
        # Start with original sequence
        optimized_seq = clean_seq
        total_changes = 0
        
        # Parse GC target
        target_gc = None
        if gc_target:
            if gc_target == 'custom':
                target_gc = 50.0  # Default custom target
            elif '-' in gc_target:
                # Range like "45-55"
                gc_range = gc_target.split('-')
                target_gc = (float(gc_range[0]) + float(gc_range[1])) / 2
        
        # Apply optimizations
        if organism and organism in self.codon_usage:
            if target_gc is not None:
                # Optimize for both GC content and organism
                optimized_seq, changes = self.optimize_gc_content(optimized_seq, target_gc, organism)
                total_changes += changes
            else:
                # Optimize for organism only
                optimized_seq, changes = self.optimize_codon_usage(optimized_seq, organism)
                total_changes += changes
        elif target_gc is not None:
            # Optimize for GC content only (use E. coli as default)
            optimized_seq, changes = self.optimize_gc_content(optimized_seq, target_gc, 'e-coli')
            total_changes += changes
        
        # Remove restriction sites if requested
        if remove_restriction:
            original_seq = optimized_seq
            optimized_seq = self.remove_restriction_sites(optimized_seq, organism or 'e-coli')
            if optimized_seq != original_seq:
                total_changes += 1
        
        # Calculate statistics
        original_gc = self.calculate_gc_content(clean_seq)
        optimized_gc = self.calculate_gc_content(optimized_seq)
        cai_score = self.calculate_cai(optimized_seq, organism or 'e-coli')
        
        # Calculate efficiency (improvement in CAI and GC target achievement)
        efficiency = 0
        if organism:
            original_cai = self.calculate_cai(clean_seq, organism)
            cai_improvement = ((cai_score - original_cai) / max(original_cai, 0.1)) * 100
            efficiency += max(0, cai_improvement)
        
        if target_gc is not None:
            original_gc_diff = abs(original_gc - target_gc)
            optimized_gc_diff = abs(optimized_gc - target_gc)
            gc_improvement = ((original_gc_diff - optimized_gc_diff) / max(original_gc_diff, 1)) * 100
            efficiency += max(0, gc_improvement)
        
        efficiency = min(100, max(0, efficiency))
        
        return {
            'success': True,
            'original_sequence': clean_seq,
            'optimized_sequence': optimized_seq,
            'protein_sequence': protein,
            'statistics': {
                'original_gc_content': round(original_gc, 2),
                'optimized_gc_content': round(optimized_gc, 2),
                'cai_score': round(cai_score, 3),
                'changes_made': total_changes,
                'efficiency_score': round(efficiency, 1),
                'sequence_length': len(optimized_seq),
                'protein_length': len(protein)
            },
            'optimization_info': {
                'gc_target': gc_target,
                'organism': organism,
                'avoid_rare_codons': avoid_rare_codons,
                'remove_restriction_sites': remove_restriction,
                'optimize_folding': optimize_folding
            }
        }

def main():
    parser = argparse.ArgumentParser(description='Codon Optimizer - Optimize DNA sequences for expression')
    parser.add_argument('--sequence', required=True, help='DNA sequence to optimize')
    parser.add_argument('--gc-target', help='GC content target (e.g., "45-55" or "custom")')
    parser.add_argument('--custom-gc', type=float, help='Custom GC target percentage (20-80)')
    parser.add_argument('--organism', help='Target organism for codon optimization')
    parser.add_argument('--avoid-rare-codons', action='store_true', default=True, help='Avoid rare codons')
    parser.add_argument('--remove-restriction-sites', action='store_true', help='Remove restriction enzyme sites')
    parser.add_argument('--optimize-folding', action='store_true', help='Optimize RNA folding')
    parser.add_argument('--file-content', help='File content for parsing')
    parser.add_argument('--file-format', help='File format (fasta or txt)')
    
    args = parser.parse_args()
    
    optimizer = CodonOptimizer()
    
    # Handle file input
    sequence = args.sequence
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
    
    # Handle custom GC target
    gc_target = args.gc_target
    if gc_target == 'custom' and args.custom_gc:
        gc_target = f"{args.custom_gc}-{args.custom_gc}"
    
    # Optimize sequence
    result = optimizer.optimize_sequence(
        sequence=sequence,
        gc_target=gc_target,
        organism=args.organism,
        avoid_rare_codons=args.avoid_rare_codons,
        remove_restriction=args.remove_restriction_sites,
        optimize_folding=args.optimize_folding
    )
    
    print(json.dumps(result))

if __name__ == '__main__':
    main()