#!/usr/bin/env python3
"""
Restriction Site Cleaner - BioPython Backend
Handles DNA sequence analysis and restriction site cleaning
"""

import sys
import json
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.MeltingTemp import Tm_NN
import re
from collections import defaultdict

# Restriction Enzyme Database
RESTRICTION_ENZYMES = {
    # 4-Base Cutters
    'MspI': 'CCGG',
    'HpaII': 'CCGG', 
    'TaqI': 'TCGA',
    'MseI': 'TTAA',
    'AluI': 'AGCT',
    'HinfI': 'GANTC',
    
    # 6-Base Cutters
    'BamHI': 'GGATCC',
    'EcoRI': 'GAATTC',
    'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG',
    'SalI': 'GTCGAC',
    'XbaI': 'TCTAGA',
    'SpeI': 'ACTAGT',
    'PstI': 'CTGCAG',
    'KpnI': 'GGTACC',
    'SacI': 'GAGCTC',
    'BglII': 'AGATCT',
    'NheI': 'GCTAGC',
    'ApaI': 'GGGCCC',
    'SmaI': 'CCCGGG',
    'EcoRV': 'GATATC',
    'HpaI': 'GTTAAC',
    
    # 8-Base Cutters
    'NotI': 'GCGGCCGC',
    'SfiI': 'GGCCNNNNNGGCC',
    'PacI': 'TTAATTAA',
    'SwaI': 'ATTTAAAT',
    
    # Type IIS Enzymes
    'BsaI': 'GGTCTC',
    'BsmBI': 'CGTCTC',
    'BbsI': 'GAAGAC',
    'Esp3I': 'CGTCTC',
    'SapI': 'GCTCTTC'
}

# Enzyme presets
ENZYME_PRESETS = {
    'common-6': ['BamHI', 'EcoRI', 'HindIII', 'XhoI', 'SalI', 'XbaI'],
    'cloning': ['BamHI', 'EcoRI', 'HindIII', 'XhoI', 'SalI', 'XbaI', 'SpeI', 'PstI', 'KpnI'],
    'blunt': ['SmaI', 'EcoRV', 'HpaI', 'AluI'],
    'type-iis': ['BsaI', 'BsmBI', 'BbsI', 'Esp3I', 'SapI']
}

# Codon usage tables for different organisms (most preferred codons for each amino acid)
CODON_TABLES = {
    'e-coli': {
        'A': ['GCG', 'GCA', 'GCC', 'GCT'],  # Alanine - GCG most preferred
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine - CGT most preferred
        'N': ['AAC', 'AAT'],  # Asparagine - AAC preferred
        'D': ['GAT', 'GAC'],  # Aspartic acid - GAT preferred
        'C': ['TGT', 'TGC'],  # Cysteine - TGT preferred
        'Q': ['CAG', 'CAA'],  # Glutamine - CAG preferred
        'E': ['GAA', 'GAG'],  # Glutamic acid - GAA preferred
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine - GGT preferred
        'H': ['CAT', 'CAC'],  # Histidine - CAT preferred
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine - ATT preferred
        'L': ['CTG', 'TTG', 'CTA', 'CTT', 'CTC', 'TTA'],  # Leucine - CTG preferred
        'K': ['AAA', 'AAG'],  # Lysine - AAA preferred
        'M': ['ATG'],  # Methionine - only ATG
        'F': ['TTT', 'TTC'],  # Phenylalanine - TTT preferred
        'P': ['CCG', 'CCA', 'CCT', 'CCC'],  # Proline - CCG preferred
        'S': ['TCG', 'AGC', 'TCA', 'TCT', 'TCC', 'AGT'],  # Serine - TCG preferred
        'T': ['ACG', 'ACA', 'ACC', 'ACT'],  # Threonine - ACG preferred
        'W': ['TGG'],  # Tryptophan - only TGG
        'Y': ['TAT', 'TAC'],  # Tyrosine - TAT preferred
        'V': ['GTG', 'GTA', 'GTT', 'GTC'],  # Valine - GTG preferred
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    },
    'human': {
        'A': ['GCC', 'GCT', 'GCA', 'GCG'],  # Alanine - GCC most preferred
        'R': ['CGC', 'CGT', 'AGA', 'AGG', 'CGA', 'CGG'],  # Arginine - CGC preferred
        'N': ['AAC', 'AAT'],  # Asparagine - AAC preferred
        'D': ['GAC', 'GAT'],  # Aspartic acid - GAC preferred
        'C': ['TGC', 'TGT'],  # Cysteine - TGC preferred
        'Q': ['CAG', 'CAA'],  # Glutamine - CAG preferred
        'E': ['GAG', 'GAA'],  # Glutamic acid - GAG preferred
        'G': ['GGC', 'GGT', 'GGA', 'GGG'],  # Glycine - GGC preferred
        'H': ['CAC', 'CAT'],  # Histidine - CAC preferred
        'I': ['ATC', 'ATT', 'ATA'],  # Isoleucine - ATC preferred
        'L': ['CTG', 'CTC', 'TTG', 'CTT', 'CTA', 'TTA'],  # Leucine - CTG preferred
        'K': ['AAG', 'AAA'],  # Lysine - AAG preferred
        'M': ['ATG'],  # Methionine - only ATG
        'F': ['TTC', 'TTT'],  # Phenylalanine - TTC preferred
        'P': ['CCC', 'CCT', 'CCA', 'CCG'],  # Proline - CCC preferred
        'S': ['TCC', 'TCT', 'AGC', 'TCA', 'TCG', 'AGT'],  # Serine - TCC preferred
        'T': ['ACC', 'ACT', 'ACA', 'ACG'],  # Threonine - ACC preferred
        'W': ['TGG'],  # Tryptophan - only TGG
        'Y': ['TAC', 'TAT'],  # Tyrosine - TAC preferred
        'V': ['GTG', 'GTC', 'GTT', 'GTA'],  # Valine - GTG preferred
        '*': ['TGA', 'TAA', 'TAG']  # Stop codons
    },
    'yeast': {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine - GCT preferred
        'R': ['AGA', 'CGT', 'CGC', 'CGA', 'CGG', 'AGG'],  # Arginine - AGA preferred
        'N': ['AAT', 'AAC'],  # Asparagine - AAT preferred
        'D': ['GAT', 'GAC'],  # Aspartic acid - GAT preferred
        'C': ['TGT', 'TGC'],  # Cysteine - TGT preferred
        'Q': ['CAA', 'CAG'],  # Glutamine - CAA preferred
        'E': ['GAA', 'GAG'],  # Glutamic acid - GAA preferred
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine - GGT preferred
        'H': ['CAT', 'CAC'],  # Histidine - CAT preferred
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine - ATT preferred
        'L': ['TTG', 'CTT', 'CTC', 'CTG', 'CTA', 'TTA'],  # Leucine - TTG preferred
        'K': ['AAA', 'AAG'],  # Lysine - AAA preferred
        'M': ['ATG'],  # Methionine - only ATG
        'F': ['TTT', 'TTC'],  # Phenylalanine - TTT preferred
        'P': ['CCT', 'CCA', 'CCC', 'CCG'],  # Proline - CCT preferred
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine - TCT preferred
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine - ACT preferred
        'W': ['TGG'],  # Tryptophan - only TGG
        'Y': ['TAT', 'TAC'],  # Tyrosine - TAT preferred
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine - GTT preferred
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    },
    'mouse': {
        # Similar to human but with slight differences
        'A': ['GCC', 'GCT', 'GCA', 'GCG'],
        'R': ['CGC', 'CGT', 'AGA', 'AGG', 'CGA', 'CGG'],
        'N': ['AAC', 'AAT'],
        'D': ['GAC', 'GAT'],
        'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'],
        'E': ['GAG', 'GAA'],
        'G': ['GGC', 'GGT', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'],
        'I': ['ATC', 'ATT', 'ATA'],
        'L': ['CTG', 'CTC', 'TTG', 'CTT', 'CTA', 'TTA'],
        'K': ['AAG', 'AAA'],
        'M': ['ATG'],
        'F': ['TTC', 'TTT'],
        'P': ['CCC', 'CCT', 'CCA', 'CCG'],
        'S': ['TCC', 'TCT', 'AGC', 'TCA', 'TCG', 'AGT'],
        'T': ['ACC', 'ACT', 'ACA', 'ACG'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
        'V': ['GTG', 'GTC', 'GTT', 'GTA'],
        '*': ['TGA', 'TAA', 'TAG']
    },
    'cho-cells': {
        # CHO cells - similar to human
        'A': ['GCC', 'GCT', 'GCA', 'GCG'],
        'R': ['CGC', 'CGT', 'AGA', 'AGG', 'CGA', 'CGG'],
        'N': ['AAC', 'AAT'],
        'D': ['GAC', 'GAT'],
        'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'],
        'E': ['GAG', 'GAA'],
        'G': ['GGC', 'GGT', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'],
        'I': ['ATC', 'ATT', 'ATA'],
        'L': ['CTG', 'CTC', 'TTG', 'CTT', 'CTA', 'TTA'],
        'K': ['AAG', 'AAA'],
        'M': ['ATG'],
        'F': ['TTC', 'TTT'],
        'P': ['CCC', 'CCT', 'CCA', 'CCG'],
        'S': ['TCC', 'TCT', 'AGC', 'TCA', 'TCG', 'AGT'],
        'T': ['ACC', 'ACT', 'ACA', 'ACG'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
        'V': ['GTG', 'GTC', 'GTT', 'GTA'],
        '*': ['TGA', 'TAA', 'TAG']
    }
}

def parse_fasta(content):
    """Parse FASTA format content"""
    lines = content.strip().split('\n')
    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    return sequence.upper()

def parse_genbank(content):
    """Parse GenBank format content"""
    lines = content.split('\n')
    sequence = ''
    in_sequence = False
    
    for line in lines:
        if line.startswith('ORIGIN'):
            in_sequence = True
            continue
        if line.startswith('//'):
            break
        if in_sequence:
            # Remove numbers and spaces from sequence lines
            clean_line = re.sub(r'[^a-zA-Z]', '', line)
            sequence += clean_line
    
    return sequence.upper()

def validate_dna_sequence(sequence):
    """Validate DNA sequence contains only valid bases"""
    # Remove whitespace
    sequence = re.sub(r'\s', '', sequence)
    # Check for valid DNA bases (including ambiguous bases)
    valid_bases = set('ATGCNRYSWKMBDHV')
    sequence_bases = set(sequence.upper())
    
    if not sequence_bases.issubset(valid_bases):
        invalid_bases = sequence_bases - valid_bases
        return False, f"Invalid bases found: {', '.join(invalid_bases)}"
    
    return True, sequence.upper()

def calculate_gc_content(sequence):
    """Calculate GC content using BioPython"""
    try:
        seq_obj = Seq(sequence)
        gc_percent = gc_fraction(seq_obj) * 100  # Convert to percentage
        return round(gc_percent, 1)
    except:
        return 0.0

def calculate_melting_temp(sequence):
    """Calculate melting temperature using BioPython"""
    try:
        seq_obj = Seq(sequence)
        # Use nearest neighbor method for accurate Tm calculation
        tm = Tm_NN(seq_obj)
        return round(tm, 1)
    except:
        return 0.0

def analyze_sequence(sequence):
    """Analyze DNA sequence and return comprehensive data"""
    try:
        # Validate sequence
        is_valid, result = validate_dna_sequence(sequence)
        if not is_valid:
            return {
                'success': False,
                'error': result
            }
        
        clean_sequence = result
        seq_obj = Seq(clean_sequence)
        
        # Basic statistics
        length = len(clean_sequence)
        gc_content = calculate_gc_content(clean_sequence)
        melting_temp = calculate_melting_temp(clean_sequence)
        
        # Base composition
        a_count = clean_sequence.count('A')
        t_count = clean_sequence.count('T')
        g_count = clean_sequence.count('G')
        c_count = clean_sequence.count('C')
        n_count = clean_sequence.count('N')
        
        # Molecular weight (approximate)
        try:
            mol_weight = molecular_weight(seq_obj, seq_type='DNA')
        except:
            mol_weight = 0
        
        return {
            'success': True,
            'data': {
                'sequence': clean_sequence,
                'length': length,
                'gc_content': gc_content,
                'melting_temp': melting_temp,
                'base_composition': {
                    'A': a_count,
                    'T': t_count,
                    'G': g_count,
                    'C': c_count,
                    'N': n_count
                },
                'molecular_weight': round(mol_weight, 2) if mol_weight else 0,
                'formatted_sequence': format_sequence_for_display(clean_sequence)
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f"Analysis error: {str(e)}"
        }

def format_sequence_for_display(sequence, line_length=140):
    """Format sequence with zero-padded line numbers like DNA Sequencer"""
    formatted = ''
    for i in range(0, len(sequence), line_length):
        line_num = i + 1
        line_seq = sequence[i:i+line_length]
        # Zero-pad line numbers to 5 digits for consistent alignment
        formatted += f"{line_num:05d}|{line_seq}\n"
    return formatted.strip()

def find_restriction_sites(sequence, enzyme_name):
    """Find all restriction sites for a given enzyme"""
    if enzyme_name not in RESTRICTION_ENZYMES:
        return []
    
    recognition_site = RESTRICTION_ENZYMES[enzyme_name]
    sites = []
    
    # Handle ambiguous bases (like N in SfiI)
    if 'N' in recognition_site:
        # For now, skip ambiguous sites - could be enhanced later
        return sites
    
    # Find all occurrences
    for i in range(len(sequence) - len(recognition_site) + 1):
        if sequence[i:i+len(recognition_site)] == recognition_site:
            sites.append({
                'position': i + 1,  # 1-based position
                'sequence': recognition_site,
                'enzyme': enzyme_name
            })
    
    return sites

def get_enzymes_to_process(enzyme_selection):
    """Get list of enzymes based on selection"""
    if enzyme_selection == 'all':
        return list(RESTRICTION_ENZYMES.keys())
    elif enzyme_selection in ENZYME_PRESETS:
        return ENZYME_PRESETS[enzyme_selection]
    elif enzyme_selection in RESTRICTION_ENZYMES:
        return [enzyme_selection]
    else:
        return []

# Standard genetic code table
GENETIC_CODE = {
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

def get_optimal_codon(amino_acid, organism, avoid_codons=None):
    """Get the most optimal codon for an amino acid in a specific organism"""
    if avoid_codons is None:
        avoid_codons = set()
    
    if organism not in CODON_TABLES:
        organism = 'e-coli'  # Default fallback
    
    codon_preferences = CODON_TABLES[organism].get(amino_acid, [])
    
    # Find the first codon that's not in the avoid list
    for codon in codon_preferences:
        if codon not in avoid_codons:
            return codon
    
    # If all preferred codons are avoided, return the first one anyway
    return codon_preferences[0] if codon_preferences else None

def optimize_codon_for_site(sequence, site_position, recognition_site, organism):
    """
    Optimize codons around a restriction site to remove it while maintaining organism preference
    """
    site_start = site_position - 1  # Convert to 0-based
    site_end = site_start + len(recognition_site)
    
    # Find the reading frame that overlaps with the restriction site
    for frame in range(3):
        # Check if this reading frame overlaps with the restriction site
        codon_start = (site_start // 3) * 3 + frame
        
        while codon_start < site_end:
            if codon_start + 3 <= len(sequence) and codon_start >= 0:
                original_codon = sequence[codon_start:codon_start + 3]
                
                # Only process if this codon overlaps with the restriction site
                if (codon_start < site_end and codon_start + 3 > site_start):
                    amino_acid = GENETIC_CODE.get(original_codon)
                    
                    if amino_acid and amino_acid != '*':  # Don't optimize stop codons
                        # Get optimal codon for this organism
                        optimal_codon = get_optimal_codon(amino_acid, organism, {original_codon})
                        
                        if optimal_codon and optimal_codon != original_codon:
                            # Replace the codon
                            new_sequence = (sequence[:codon_start] + 
                                          optimal_codon + 
                                          sequence[codon_start + 3:])
                            
                            # Check if this removes the restriction site
                            if recognition_site not in new_sequence[max(0, site_start-10):site_end+10]:
                                return new_sequence, {
                                    'position': codon_start + 1,
                                    'original_codon': original_codon,
                                    'new_codon': optimal_codon,
                                    'amino_acid': amino_acid,
                                    'change_type': 'codon_optimization'
                                }
            
            codon_start += 3
    
    return sequence, None

def create_silent_mutation(codon, target_aa):
    """Create a silent mutation for a codon to avoid restriction site"""
    # Find alternative codons for the same amino acid
    alternatives = [c for c, aa in GENETIC_CODE.items() if aa == target_aa and c != codon]
    return alternatives[0] if alternatives else codon

def clean_restriction_sites(sequence, enzymes, options=None):
    """Clean restriction sites from sequence"""
    if options is None:
        options = {
            'preserve_orfs': True,
            'silent_mutations': True,
            'optimize_codons': False,
            'preserve_regulatory': True
        }
    
    cleaned_sequence = sequence
    changes_made = []
    sites_found = defaultdict(list)
    
    # Find all restriction sites
    for enzyme in enzymes:
        sites = find_restriction_sites(sequence, enzyme)
        if sites:
            sites_found[enzyme] = sites
    
    # Process each enzyme's sites
    total_sites_removed = 0
    for enzyme, sites in sites_found.items():
        recognition_site = RESTRICTION_ENZYMES[enzyme]
        
        for site in sites:
            pos = site['position'] - 1  # Convert to 0-based
            
            # Check if this site still exists in the current cleaned sequence
            if pos < len(cleaned_sequence) - len(recognition_site) + 1:
                current_site = cleaned_sequence[pos:pos+len(recognition_site)]
                
                if current_site == recognition_site:
                    mutation_applied = False
                    
                    # Try codon optimization first if enabled
                    if options.get('optimize_codons') and options.get('target_organism'):
                        new_sequence, codon_change = optimize_codon_for_site(
                            cleaned_sequence, site['position'], recognition_site, options['target_organism']
                        )
                        
                        if codon_change:
                            cleaned_sequence = new_sequence
                            changes_made.append({
                                'enzyme': enzyme,
                                'position': site['position'],
                                'original': codon_change['original_codon'],
                                'mutated': codon_change['new_codon'],
                                'amino_acid': codon_change['amino_acid'],
                                'change_type': 'codon_optimization',
                                'organism': options['target_organism']
                            })
                            total_sites_removed += 1
                            mutation_applied = True
                    
                    # If codon optimization didn't work, try simple mutation
                    if not mutation_applied:
                        # Simple mutation strategy: change one base
                        mutated_site = list(current_site)
                        
                        # Change the first base to something different
                        bases = ['A', 'T', 'G', 'C']
                        for base in bases:
                            if base != mutated_site[0]:
                                mutated_site[0] = base
                                break
                        
                        new_site = ''.join(mutated_site)
                        
                        # Apply the mutation
                        cleaned_sequence = (cleaned_sequence[:pos] + 
                                          new_site + 
                                          cleaned_sequence[pos+len(recognition_site):])
                        
                        changes_made.append({
                            'enzyme': enzyme,
                            'position': site['position'],
                            'original': current_site,
                            'mutated': new_site,
                            'change_type': 'silent_mutation' if options.get('silent_mutations') else 'base_substitution'
                        })
                        
                        total_sites_removed += 1
    
    return {
        'cleaned_sequence': cleaned_sequence,
        'original_sequence': sequence,
        'sites_found': dict(sites_found),
        'changes_made': changes_made,
        'statistics': {
            'total_sites_found': sum(len(sites) for sites in sites_found.values()),
            'sites_removed': total_sites_removed,
            'enzymes_processed': len(sites_found),
            'sequence_length_original': len(sequence),
            'sequence_length_cleaned': len(cleaned_sequence)
        }
    }

def comprehensive_site_scan(sequence):
    """
    Comprehensive scan for ALL restriction sites in the sequence
    Returns detailed information about every site found
    """
    sites_found = {}
    total_sites = 0
    enzymes_found = 0
    
    # Scan for all enzymes in the database
    for enzyme_name in RESTRICTION_ENZYMES.keys():
        sites = find_restriction_sites(sequence, enzyme_name)
        
        if sites:
            sites_found[enzyme_name] = sites
            total_sites += len(sites)
            enzymes_found += 1
    
    return {
        'sites_found': sites_found,
        'total_sites': total_sites,
        'enzymes_found': enzymes_found
    }

def targeted_site_scan(sequence, enzymes_to_scan):
    """
    Targeted scan for specific restriction enzymes in the sequence
    Returns detailed information about sites found for the specified enzymes
    """
    sites_found = {}
    total_sites = 0
    enzymes_found = 0
    
    # Scan only for specified enzymes
    for enzyme_name in enzymes_to_scan:
        if enzyme_name in RESTRICTION_ENZYMES:
            sites = find_restriction_sites(sequence, enzyme_name)
            
            if sites:
                sites_found[enzyme_name] = sites
                total_sites += len(sites)
                enzymes_found += 1
    
    return {
        'sites_found': sites_found,
        'total_sites': total_sites,
        'enzymes_found': enzymes_found
    }

def main():
    """Main function to handle command line arguments"""
    try:
        if len(sys.argv) < 2:
            print(json.dumps({
                'success': False,
                'error': 'No operation specified'
            }))
            return
        
        operation = sys.argv[1]
        
        if operation == 'load':
            # Load and analyze sequence
            if len(sys.argv) < 3:
                print(json.dumps({
                    'success': False,
                    'error': 'No sequence provided'
                }))
                return
            
            sequence_input = sys.argv[2]
            
            # Check if it's a file format or raw sequence
            if sequence_input.startswith('>'):
                # FASTA format
                sequence = parse_fasta(sequence_input)
            elif 'ORIGIN' in sequence_input:
                # GenBank format
                sequence = parse_genbank(sequence_input)
            else:
                # Raw sequence
                sequence = sequence_input.strip()
            
            # Analyze the sequence
            result = analyze_sequence(sequence)
            print(json.dumps(result))
            
        elif operation == 'analyze':
            # Analyze and clean restriction sites
            if len(sys.argv) < 4:
                print(json.dumps({
                    'success': False,
                    'error': 'Insufficient arguments for analysis'
                }))
                return
            
            sequence_input = sys.argv[2]
            enzyme_selection = sys.argv[3]
            
            # Parse options if provided
            options = {
                'preserve_orfs': True,
                'silent_mutations': True,
                'optimize_codons': False,
                'preserve_regulatory': True
            }
            
            if len(sys.argv) > 4:
                try:
                    options.update(json.loads(sys.argv[4]))
                except:
                    pass  # Use defaults if parsing fails
            
            # Parse sequence
            if sequence_input.startswith('>'):
                sequence = parse_fasta(sequence_input)
            elif 'ORIGIN' in sequence_input:
                sequence = parse_genbank(sequence_input)
            else:
                sequence = sequence_input.strip()
            
            # Validate sequence
            is_valid, result = validate_dna_sequence(sequence)
            if not is_valid:
                print(json.dumps({
                    'success': False,
                    'error': result
                }))
                return
            
            clean_sequence = result
            
            # Get enzymes to process
            enzymes = get_enzymes_to_process(enzyme_selection)
            if not enzymes:
                print(json.dumps({
                    'success': False,
                    'error': f'Unknown enzyme selection: {enzyme_selection}'
                }))
                return
            
            # Clean restriction sites
            cleaning_result = clean_restriction_sites(clean_sequence, enzymes, options)
            
            # Analyze cleaned sequence
            cleaned_analysis = analyze_sequence(cleaning_result['cleaned_sequence'])
            
            if cleaned_analysis['success']:
                # Combine results
                final_result = {
                    'success': True,
                    'cleaned_sequence': cleaning_result['cleaned_sequence'],
                    'original_sequence': cleaning_result['original_sequence'],
                    'sites_found': cleaning_result['sites_found'],
                    'changes_made': cleaning_result['changes_made'],
                    'statistics': cleaning_result['statistics'],
                    'sequence_analysis': cleaned_analysis['data'],
                    'enzymes_processed': enzymes,
                    'options_used': options
                }
                
                print(json.dumps(final_result))
            else:
                print(json.dumps({
                    'success': False,
                    'error': f'Failed to analyze cleaned sequence: {cleaned_analysis.get("error", "Unknown error")}'
                }))
            
        elif operation == 'scan':
            # Comprehensive restriction site analysis (no cleaning)
            if len(sys.argv) < 3:
                print(json.dumps({
                    'success': False,
                    'error': 'No sequence provided for scanning'
                }))
                return
            
            sequence_input = sys.argv[2]
            
            # Optional enzyme selection for targeted scanning
            enzyme_selection = None
            if len(sys.argv) > 3:
                enzyme_selection = sys.argv[3]
            
            # Parse sequence
            if sequence_input.startswith('>'):
                sequence = parse_fasta(sequence_input)
            elif 'ORIGIN' in sequence_input:
                sequence = parse_genbank(sequence_input)
            else:
                sequence = sequence_input.strip()
            
            # Validate sequence
            is_valid, result = validate_dna_sequence(sequence)
            if not is_valid:
                print(json.dumps({
                    'success': False,
                    'error': result
                }))
                return
            
            clean_sequence = result
            
            # Scan for restriction sites (all or specific enzymes)
            if enzyme_selection:
                # Get specific enzymes to scan
                enzymes_to_scan = get_enzymes_to_process(enzyme_selection)
                if not enzymes_to_scan:
                    print(json.dumps({
                        'success': False,
                        'error': f'Unknown enzyme selection: {enzyme_selection}'
                    }))
                    return
                scan_result = targeted_site_scan(clean_sequence, enzymes_to_scan)
            else:
                # Scan all enzymes
                scan_result = comprehensive_site_scan(clean_sequence)
            
            # Analyze original sequence
            sequence_analysis = analyze_sequence(clean_sequence)
            
            if sequence_analysis['success']:
                # Combine results
                final_result = {
                    'success': True,
                    'sequence': clean_sequence,
                    'sites_found': scan_result['sites_found'],
                    'total_sites': scan_result['total_sites'],
                    'enzymes_found': scan_result['enzymes_found'],
                    'sequence_analysis': sequence_analysis['data']
                }
                
                print(json.dumps(final_result))
            else:
                print(json.dumps({
                    'success': False,
                    'error': f'Failed to analyze sequence: {sequence_analysis.get("error", "Unknown error")}'
                }))
            
        else:
            print(json.dumps({
                'success': False,
                'error': f'Unknown operation: {operation}'
            }))
            
    except Exception as e:
        print(json.dumps({
            'success': False,
            'error': f'Script error: {str(e)}'
        }))

if __name__ == '__main__':
    main()