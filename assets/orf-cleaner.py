#!/usr/bin/env python3
"""
ORF Cleaner - BioPython Backend
Handles DNA sequence analysis and ORF (Open Reading Frame) optimization
"""

import sys
import json
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.MeltingTemp import Tm_NN
import re
from collections import defaultdict

# Standard genetic code
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
        'L': ['CTG', 'CTC', 'TTG', 'TTA', 'CTT', 'CTA'],  # Leucine - CTG preferred
        'K': ['AAG', 'AAA'],  # Lysine - AAG preferred
        'M': ['ATG'],  # Methionine - only ATG
        'F': ['TTC', 'TTT'],  # Phenylalanine - TTC preferred
        'P': ['CCC', 'CCT', 'CCA', 'CCG'],  # Proline - CCC preferred
        'S': ['AGC', 'TCC', 'TCT', 'TCA', 'TCG', 'AGT'],  # Serine - AGC preferred
        'T': ['ACC', 'ACT', 'ACA', 'ACG'],  # Threonine - ACC preferred
        'W': ['TGG'],  # Tryptophan - only TGG
        'Y': ['TAC', 'TAT'],  # Tyrosine - TAC preferred
        'V': ['GTG', 'GTC', 'GTT', 'GTA'],  # Valine - GTG preferred
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    },
    'yeast': {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine - GCT most preferred
        'R': ['AGA', 'CGT', 'CGC', 'CGA', 'CGG', 'AGG'],  # Arginine - AGA preferred
        'N': ['AAT', 'AAC'],  # Asparagine - AAT preferred
        'D': ['GAT', 'GAC'],  # Aspartic acid - GAT preferred
        'C': ['TGT', 'TGC'],  # Cysteine - TGT preferred
        'Q': ['CAA', 'CAG'],  # Glutamine - CAA preferred
        'E': ['GAA', 'GAG'],  # Glutamic acid - GAA preferred
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine - GGT preferred
        'H': ['CAT', 'CAC'],  # Histidine - CAT preferred
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine - ATT preferred
        'L': ['TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA'],  # Leucine - TTG preferred
        'K': ['AAA', 'AAG'],  # Lysine - AAA preferred
        'M': ['ATG'],  # Methionine - only ATG
        'F': ['TTT', 'TTC'],  # Phenylalanine - TTT preferred
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline - CCT preferred
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine - TCT preferred
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine - ACT preferred
        'W': ['TGG'],  # Tryptophan - only TGG
        'Y': ['TAT', 'TAC'],  # Tyrosine - TAT preferred
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine - GTT preferred
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
}

def validate_dna_sequence(sequence):
    """Validate and clean DNA sequence"""
    if not sequence:
        return False, "Empty sequence provided"
    
    # Remove whitespace and convert to uppercase
    clean_seq = re.sub(r'\s+', '', sequence.upper())
    
    # Check for invalid characters
    invalid_chars = re.findall(r'[^ATGCNRYSWKMBDHV]', clean_seq)
    if invalid_chars:
        return False, f"Invalid DNA characters found: {', '.join(set(invalid_chars))}"
    
    # Convert ambiguous bases to N for processing
    clean_seq = re.sub(r'[RYSWKMBDHV]', 'N', clean_seq)
    
    if len(clean_seq) < 3:
        return False, "Sequence too short (minimum 3 nucleotides required)"
    
    return True, clean_seq

def find_orfs(sequence, min_length=30):
    """Find Open Reading Frames - SIMPLE SINGLE FRAME MODE
    
    1. Find first ATG (start codon)
    2. Read 3 bases at a time from there
    3. Flag any additional ATG found after the first one
    4. Flag any stop codon found in the middle (not the final one)
    5. Continue until end or stop codon
    """
    orfs = []
    
    # Make sure sequence length is divisible by 3
    clean_seq = sequence
    while len(clean_seq) % 3 != 0:
        clean_seq = clean_seq[:-1]
    
    # Find the FIRST ATG in the sequence
    first_atg_pos = -1
    for i in range(0, len(clean_seq) - 2, 3):
        codon = clean_seq[i:i+3]
        if codon == 'ATG':
            first_atg_pos = i
            break
    
    # If no ATG found, return empty
    if first_atg_pos == -1:
        return orfs
    
    # Start reading from the first ATG
    start_pos = first_atg_pos
    end_pos = len(clean_seq)  # Read to the very end of sequence
    
    # Extract the ENTIRE ORF sequence from first ATG to end
    orf_seq = clean_seq[start_pos:end_pos]
    orf_length = end_pos - start_pos
    protein = translate_sequence(orf_seq)
    
    # Create the single ORF
    orfs.append({
        'start': start_pos + 1,  # 1-based indexing
        'end': end_pos,
        'length': orf_length,
        'frame': 1,
        'sequence': orf_seq,
        'protein': protein,
        'protein_length': len(protein) - 1 if protein.endswith('*') else len(protein),
        'meets_min_length': orf_length >= min_length
    })
    
    return orfs

def translate_sequence(dna_sequence):
    """Translate DNA sequence to protein"""
    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = GENETIC_CODE.get(codon, 'X')
            protein += amino_acid
    return protein

def fix_problematic_codons(sequence, organism='e-coli', fix_internal_stops=True, fix_internal_starts=True):
    """Fix problematic codons (internal stops and starts) by mutating them to organism-preferred codons"""
    if organism not in CODON_TABLES:
        return sequence, 0, []
    
    codon_table = CODON_TABLES[organism]
    fixed_seq = ""
    changes_made = 0
    fix_log = []
    
    stop_codons = ['TAA', 'TAG', 'TGA']
    start_codons = ['ATG']
    
    # Process sequence in codons
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            should_fix = False
            fix_reason = ""
            
            # Check if this is an internal stop codon (not the last codon)
            if fix_internal_stops and codon in stop_codons and i < len(sequence) - 3:
                should_fix = True
                fix_reason = "internal_stop"
            
            # Check if this is an internal start codon (not the first codon)
            elif fix_internal_starts and codon in start_codons and i > 0:
                should_fix = True
                fix_reason = "internal_start"
            
            if should_fix:
                # Find a suitable replacement codon
                # Priority: Leucine (L) -> Serine (S) -> Alanine (A) -> Glycine (G)
                replacement_amino_acids = ['L', 'S', 'A', 'G']
                replacement_codon = None
                replacement_aa = None
                
                for aa in replacement_amino_acids:
                    if aa in codon_table and codon_table[aa]:
                        replacement_codon = codon_table[aa][0]  # Most preferred codon
                        replacement_aa = aa
                        break
                
                if replacement_codon:
                    fix_log.append({
                        'position': i + 1,
                        'original_codon': codon,
                        'fixed_codon': replacement_codon,
                        'original_amino_acid': GENETIC_CODE.get(codon, 'X'),
                        'new_amino_acid': replacement_aa,
                        'reason': fix_reason
                    })
                    changes_made += 1
                    fixed_seq += replacement_codon
                else:
                    # Fallback: use the original codon if no replacement found
                    fixed_seq += codon
            else:
                fixed_seq += codon
        else:
            fixed_seq += codon
    
    return fixed_seq, changes_made, fix_log

def optimize_codon_usage(sequence, organism='e-coli'):
    """Smart codon optimization - only fix problematic codons, keep good ones"""
    if organism not in CODON_TABLES:
        return sequence, 0, []
    
    codon_table = CODON_TABLES[organism]
    optimized_seq = ""
    changes_made = 0
    optimization_log = []
    
    # Define problematic codons that should be replaced
    STOP_CODONS = {'TGA', 'TAA', 'TAG'}
    
    # Define rare/problematic codons for each organism that should be avoided
    # These are codons that are significantly less preferred and should be replaced
    RARE_CODONS = {
        'e-coli': {
            'CGA', 'CGG', 'AGA', 'AGG',  # Rare arginine codons (prefer CGT, CGC)
            'CTA', 'TTA',  # Rare leucine codons (prefer CTG)
            'ATA',  # Rare isoleucine codon (prefer ATT, ATC)
            'GGA', 'GGG',  # Less preferred glycine (prefer GGT, GGC)
            'CCC', 'CCA', 'CCT',  # Less preferred proline (prefer CCG)
            'TCA', 'AGT',  # Less preferred serine (prefer TCG, AGC)
            'GCA', 'GCT', 'GCC',  # Less preferred alanine (prefer GCG)
        },
        'yeast': {
            'CGG', 'CGA', 'CGC', 'AGG',  # Rare arginine (prefer AGA, CGT)
            'CCG', 'CCA', 'CCC',  # Rare proline (prefer CCT)
            'GCG', 'GCA', 'GCC',  # Rare alanine (prefer GCT)
            'TCG', 'TCA', 'AGC',  # Rare serine (prefer TCT, TCC)
            'CTG', 'CTA', 'TTA',  # Less preferred leucine (prefer TTG)
            'ACG', 'ACA', 'ACC',  # Less preferred threonine (prefer ACT)
        },
        'human': {
            'CGA', 'CGG',  # Rare arginine codons (prefer CGC, CGT)
            'CCG',  # Rare proline (prefer CCC, CCT)
            'GCG',  # Rare alanine (prefer GCC, GCT)
            'TCG',  # Rare serine (prefer AGC, TCC)
            'TTA', 'CTA',  # Less preferred leucine (prefer CTG, CTC)
            'ACG',  # Less preferred threonine (prefer ACC, ACT)
            'GTA',  # Less preferred valine (prefer GTG, GTC)
        }
    }
    
    rare_codons_for_organism = RARE_CODONS.get(organism, set())
    
    # Process sequence in codons
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = GENETIC_CODE.get(codon, 'X')
            
            # Only optimize if codon is problematic
            should_optimize = False
            reason = ""
            
            if codon in STOP_CODONS:
                should_optimize = True
                reason = "stop_codon"
            elif codon in rare_codons_for_organism:
                should_optimize = True
                reason = "rare_codon"
            
            if should_optimize and amino_acid != 'X' and amino_acid in codon_table:
                # Get the most preferred codon for this amino acid
                preferred_codon = codon_table[amino_acid][0]
                
                optimization_log.append({
                    'position': i + 1,
                    'original_codon': codon,
                    'optimized_codon': preferred_codon,
                    'amino_acid': amino_acid,
                    'reason': reason
                })
                changes_made += 1
                optimized_seq += preferred_codon
            else:
                # Keep the original codon (it's fine)
                optimized_seq += codon
        else:
            optimized_seq += codon
    
    return optimized_seq, changes_made, optimization_log

def analyze_sequence(sequence):
    """Analyze DNA sequence properties"""
    try:
        # Validate sequence
        is_valid, clean_seq = validate_dna_sequence(sequence)
        if not is_valid:
            return {
                'success': False,
                'error': clean_seq
            }
        
        # Basic sequence properties
        length = len(clean_seq)
        gc_content = round(gc_fraction(clean_seq) * 100, 1)
        
        # Calculate melting temperature
        try:
            melting_temp = round(Tm_NN(Seq(clean_seq)), 1)
        except:
            melting_temp = 0.0
        
        # Format sequence for display (with line numbers)
        formatted_sequence = format_sequence_for_display(clean_seq)
        
        return {
            'success': True,
            'data': {
                'sequence': clean_seq,
                'formatted_sequence': formatted_sequence,
                'length': length,
                'gc_content': gc_content,
                'melting_temp': melting_temp
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f'Analysis failed: {str(e)}'
        }

def preview_orf_issues(sequence, options):
    """Preview ORF issues without fixing them"""
    try:
        # Find ORFs
        orfs = find_orfs(sequence, options.get('min_orf_length', 30))
        
        issues_found = []
        total_issues = 0
        
        for orf in orfs:
            orf_sequence = orf['sequence']
            orf_issues = []
            
            # Check for internal/premature stop codons
            if options.get('remove_internal_stops', True):
                stop_codons = ['TAA', 'TAG', 'TGA']
                # Check all codons after start codon, but exclude the final stop codon
                for i in range(3, len(orf_sequence) - 3, 3):  # Skip start, check until before final stop
                    if i + 3 <= len(orf_sequence):  # Make sure we have a complete codon
                        codon = orf_sequence[i:i+3]
                        if len(codon) == 3 and codon in stop_codons:
                            position_in_orf = i + 1  # 1-based position within ORF
                            global_position = orf['start'] + i
                            
                            if position_in_orf <= 100:
                                # Premature stop codon (within first 100 bp)
                                orf_issues.append({
                                    'type': 'premature_stop_codon',
                                    'position': global_position,
                                    'orf_position': position_in_orf,
                                    'codon': codon,
                                    'description': f'Premature stop codon {codon} at position {position_in_orf} bp in ORF (global position {global_position})'
                                })
                            else:
                                # Regular internal stop codon (after 100 bp)
                                orf_issues.append({
                                    'type': 'internal_stop_codon',
                                    'position': global_position,
                                    'orf_position': position_in_orf,
                                    'codon': codon,
                                    'description': f'Internal stop codon {codon} at position {position_in_orf} bp in ORF (global position {global_position})'
                                })
                
                # Check for very short ORFs (premature termination)
                if len(orf_sequence) <= 30:  # ORFs shorter than 30 bp (10 amino acids) are likely premature
                    final_codon = orf_sequence[-3:]  # Get the final stop codon
                    if final_codon in stop_codons:
                        orf_issues.append({
                            'type': 'premature_stop_codon',
                            'position': orf['start'] + len(orf_sequence) - 3,
                            'orf_position': len(orf_sequence) - 2,  # Position within ORF (1-based)
                            'codon': final_codon,
                            'description': f'Premature stop codon {final_codon} - ORF too short ({len(orf_sequence)} bp, {len(orf_sequence)//3} amino acids) at global position {orf["start"] + len(orf_sequence) - 3}'
                        })
            
            # Check for internal start codons (ATG in middle of ORF)
            start_codons = ['ATG']
            for i in range(3, len(orf_sequence) - 3, 3):  # Skip the first start codon
                codon = orf_sequence[i:i+3]
                if codon in start_codons:
                    position_in_orf = i + 1  # 1-based position within ORF
                    global_position = orf['start'] + i
                    orf_issues.append({
                        'type': 'internal_start_codon',
                        'position': global_position,
                        'orf_position': position_in_orf,
                        'codon': codon,
                        'description': f'Internal start codon {codon} at position {position_in_orf} bp in ORF (global position {global_position})'
                    })
            
            if orf_issues:
                issues_found.append({
                    'orf_id': len(issues_found) + 1,
                    'orf_start': orf['start'],
                    'orf_end': orf['end'],
                    'orf_length': orf['length'],
                    'issues': orf_issues,
                    'issues_count': len(orf_issues)
                })
                total_issues += len(orf_issues)
        
        # Analyze original sequence
        sequence_analysis = analyze_sequence(sequence)
        
        return {
            'success': True,
            'preview_mode': True,
            'original_sequence': sequence,
            'formatted_sequence': format_sequence_for_display(sequence),
            'orfs_found': len(orfs),
            'orfs_with_issues': len(issues_found),
            'total_issues': total_issues,
            'issues_found': issues_found,
            'sequence_analysis': sequence_analysis['data'] if sequence_analysis['success'] else None,
            'options_used': options,
            'statistics': {
                'orfs_found': len(orfs),
                'orfs_with_issues': len(issues_found),
                'total_issues': total_issues,
                'issues_by_type': {
                    'premature_stop_codons': sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == 'premature_stop_codon'),
                    'internal_stop_codons': sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == 'internal_stop_codon'),
                    'internal_start_codons': sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == 'internal_start_codon'),
                    'suboptimal_codons': sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == 'suboptimal_codon')
                },
                'detailed_issues': [
                    {
                        'type': issue_type,
                        'count': sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == issue_type),
                        'locations': [
                            {
                                'orf_id': orf['orf_id'],
                                'position': issue['position'],
                                'orf_position': issue.get('orf_position', issue['position'] - orf['orf_start'] + 1),
                                'codon': issue.get('codon', issue.get('current_codon', '')),
                                'description': issue['description']
                            }
                            for orf in issues_found 
                            for issue in orf['issues'] 
                            if issue['type'] == issue_type
                        ]
                    }
                    for issue_type in ['premature_stop_codon', 'internal_stop_codon', 'internal_start_codon', 'suboptimal_codon']
                    if sum(1 for orf in issues_found for issue in orf['issues'] if issue['type'] == issue_type) > 0
                ]
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f'Preview failed: {str(e)}'
        }

def format_sequence_for_display(sequence, line_length=130):
    """Format sequence with line numbers for display"""
    formatted_lines = []
    for i in range(0, len(sequence), line_length):
        line_num = str(i + 1).zfill(5)
        line_seq = sequence[i:i + line_length]
        formatted_lines.append(f"{line_num}| {line_seq}")
    return '\n'.join(formatted_lines)

def clean_orfs(sequence, options=None):
    """Find and optimize ORFs in the sequence"""
    if options is None:
        options = {
            'min_orf_length': 30,
            'optimize_codons': True,
            'target_organism': 'e-coli',
            'preserve_start_codons': True,
            'remove_internal_stops': True
        }
    
    try:
        # Validate sequence
        is_valid, clean_seq = validate_dna_sequence(sequence)
        if not is_valid:
            return {
                'success': False,
                'error': clean_seq
            }
        
        # Find ORFs
        min_length = options.get('min_orf_length', 30)
        orfs_found = find_orfs(clean_seq, min_length)
        
        if not orfs_found:
            return {
                'success': True,
                'cleaned_sequence': clean_seq,
                'original_sequence': clean_seq,
                'orfs_found': [],
                'orfs_optimized': 0,
                'codons_improved': 0,
                'optimization_log': [],
                'statistics': {
                    'orfs_found': 0,
                    'orfs_optimized': 0,
                    'codons_improved': 0,
                    'total_improvements': 0
                }
            }
        
        # Optimize ORFs based on selected options
        optimized_sequence = clean_seq
        total_codons_improved = 0
        total_fixes_made = 0
        optimization_log = []
        fix_log = []
        orfs_optimized = 0
        
        organism = options.get('target_organism', 'e-coli')
        
        # Process each ORF
        for orf in orfs_found:
            orf_seq = orf['sequence']
            orf_modified = False
            current_orf_seq = orf_seq
            
            # Step 1: Fix problematic codons (internal stops/starts) if requested
            fix_internal_stops = options.get('remove_internal_stops', True)
            fix_internal_starts = options.get('fix_internal_starts', True)
            
            if (fix_internal_stops or fix_internal_starts) and organism:
                fixed_orf, fixes_made, orf_fix_log = fix_problematic_codons(
                    current_orf_seq, organism, fix_internal_stops, fix_internal_starts
                )
                if fixes_made > 0:
                    current_orf_seq = fixed_orf
                    total_fixes_made += fixes_made
                    fix_log.extend(orf_fix_log)
                    orf_modified = True
            
            # Step 2: Optimize codon usage if requested
            if options.get('optimize_codons', True) and organism:
                optimized_orf, codons_improved, orf_opt_log = optimize_codon_usage(current_orf_seq, organism)
                if codons_improved > 0:
                    current_orf_seq = optimized_orf
                    total_codons_improved += codons_improved
                    optimization_log.extend(orf_opt_log)
                    orf_modified = True
            
            # Replace the ORF in the main sequence if any modifications were made
            if orf_modified:
                start_idx = orf['start'] - 1  # Convert to 0-based
                end_idx = orf['end']
                optimized_sequence = (optimized_sequence[:start_idx] + 
                                    current_orf_seq + 
                                    optimized_sequence[end_idx:])
                orfs_optimized += 1
        
        return {
            'success': True,
            'cleaned_sequence': optimized_sequence,
            'original_sequence': clean_seq,
            'orfs_found': orfs_found,
            'orfs_optimized': orfs_optimized,
            'codons_improved': total_codons_improved,
            'codons_fixed': total_fixes_made,
            'optimization_log': optimization_log,
            'fix_log': fix_log,
            'statistics': {
                'orfs_found': len(orfs_found),
                'orfs_optimized': orfs_optimized,
                'codons_improved': total_codons_improved,
                'codons_fixed': total_fixes_made,
                'total_improvements': orfs_optimized + total_codons_improved + total_fixes_made
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f'ORF cleaning failed: {str(e)}'
        }

def parse_fasta(content):
    """Parse FASTA format content"""
    lines = content.split('\n')
    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    return sequence

def parse_genbank(content):
    """Parse GenBank format content"""
    lines = content.split('\n')
    sequence = ''
    in_sequence = False
    
    for line in lines:
        if line.startswith('ORIGIN'):
            in_sequence = True
            continue
        elif line.startswith('//'):
            break
        elif in_sequence:
            # Remove line numbers and spaces
            seq_line = re.sub(r'^\s*\d+', '', line)
            seq_line = re.sub(r'\s+', '', seq_line)
            sequence += seq_line
    
    return sequence

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
            # Analyze ORFs (preview mode - find issues but don't fix them)
            if len(sys.argv) < 3:
                print(json.dumps({
                    'success': False,
                    'error': 'No sequence provided for ORF analysis'
                }))
                return
            
            sequence_input = sys.argv[2]
            
            # Parse options if provided
            options = {
                'min_orf_length': 30,
                'optimize_codons': True,
                'target_organism': 'e-coli',
                'preserve_start_codons': True,
                'remove_internal_stops': True
            }
            
            if len(sys.argv) > 3:
                try:
                    parsed_options = json.loads(sys.argv[3])
                    options.update(parsed_options)
                except Exception as e:
                    pass  # Use defaults if parsing fails
            
            # Parse sequence
            if sequence_input.startswith('>'):
                sequence = parse_fasta(sequence_input)
            elif 'ORIGIN' in sequence_input:
                sequence = parse_genbank(sequence_input)
            else:
                sequence = sequence_input.strip()
            
            # Preview ORF issues (don't actually fix them)
            preview_result = preview_orf_issues(sequence, options)
            print(json.dumps(preview_result))
            
        elif operation == 'optimize':
            # Actually optimize and clean ORFs
            if len(sys.argv) < 3:
                print(json.dumps({
                    'success': False,
                    'error': 'No sequence provided for ORF optimization'
                }))
                return
            
            sequence_input = sys.argv[2]
            
            # Parse options if provided
            options = {
                'min_orf_length': 30,
                'optimize_codons': True,
                'target_organism': 'e-coli',
                'preserve_start_codons': True,
                'remove_internal_stops': True
            }
            
            if len(sys.argv) > 3:
                try:
                    options.update(json.loads(sys.argv[3]))
                except:
                    pass  # Use defaults if parsing fails
            
            # Parse sequence
            if sequence_input.startswith('>'):
                sequence = parse_fasta(sequence_input)
            elif 'ORIGIN' in sequence_input:
                sequence = parse_genbank(sequence_input)
            else:
                sequence = sequence_input.strip()
            
            # Clean ORFs
            cleaning_result = clean_orfs(sequence, options)
            
            if cleaning_result['success']:
                # Analyze cleaned sequence
                cleaned_analysis = analyze_sequence(cleaning_result['cleaned_sequence'])
                
                if cleaned_analysis['success']:
                    # Combine results
                    final_result = {
                        'success': True,
                        'cleaned_sequence': cleaning_result['cleaned_sequence'],
                        'original_sequence': cleaning_result['original_sequence'],
                        'orfs_found': cleaning_result['orfs_found'],
                        'orfs_optimized': cleaning_result['orfs_optimized'],
                        'codons_improved': cleaning_result['codons_improved'],
                        'optimization_log': cleaning_result['optimization_log'],
                        'statistics': cleaning_result['statistics'],
                        'sequence_analysis': cleaned_analysis['data'],
                        'options_used': options
                    }
                    
                    print(json.dumps(final_result))
                else:
                    print(json.dumps({
                        'success': False,
                        'error': f'Failed to analyze cleaned sequence: {cleaned_analysis.get("error", "Unknown error")}'
                    }))
            else:
                print(json.dumps(cleaning_result))
            
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