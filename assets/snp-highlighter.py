#!/usr/bin/env python3
"""
SNP Highlighter - BioPython Backend
Handles DNA sequence analysis and SNP detection
"""

import sys
import json
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.MeltingTemp import Tm_NN
import re
from collections import defaultdict

def validate_dna_sequence(sequence):
    """Validate and clean DNA sequence"""
    if not sequence:
        return False, "Empty sequence provided"
    
    # Remove whitespace and convert to uppercase
    clean_sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Check for valid DNA bases (including ambiguous bases)
    valid_bases = set('ATGCNRYSWKMBDHV')
    invalid_bases = set(clean_sequence) - valid_bases
    
    if invalid_bases:
        return False, f"Invalid DNA bases found: {', '.join(sorted(invalid_bases))}"
    
    if len(clean_sequence) == 0:
        return False, "No valid DNA sequence found"
    
    return True, clean_sequence

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    if not sequence:
        return 0.0
    
    try:
        seq_obj = Seq(sequence)
        gc_content = gc_fraction(seq_obj) * 100
        return round(gc_content, 1)
    except:
        # Fallback calculation
        gc_count = sequence.count('G') + sequence.count('C')
        return round((gc_count / len(sequence)) * 100, 1) if len(sequence) > 0 else 0.0

def calculate_melting_temp(sequence):
    """Calculate melting temperature"""
    if not sequence or len(sequence) < 10:
        return 0.0
    
    try:
        seq_obj = Seq(sequence)
        # Use nearest neighbor method for more accurate results
        tm = Tm_NN(seq_obj)
        return round(tm, 1)
    except:
        # Fallback: simple calculation
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        tm = (gc_count * 4) + (at_count * 2)
        return round(tm, 1)

def analyze_sequence(sequence, seq_type='unknown'):
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
                'type': seq_type,
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

def format_sequence_for_display(sequence, line_length=80):
    """Format sequence with zero-padded line numbers like DNA Sequencer"""
    formatted = ''
    for i in range(0, len(sequence), line_length):
        line_num = i + 1
        line_seq = sequence[i:i+line_length]
        # Zero-pad line numbers to 5 digits for consistent alignment
        formatted += f"{line_num:05d}|{line_seq}\n"
    return formatted.strip()

def find_snps(reference_seq, sample_seq):
    """Find Single Nucleotide Polymorphisms between two sequences"""
    try:
        # Validate both sequences
        is_valid_ref, ref_result = validate_dna_sequence(reference_seq)
        if not is_valid_ref:
            return {
                'success': False,
                'error': f"Reference sequence error: {ref_result}"
            }
        
        is_valid_sample, sample_result = validate_dna_sequence(sample_seq)
        if not is_valid_sample:
            return {
                'success': False,
                'error': f"Sample sequence error: {sample_result}"
            }
        
        ref_clean = ref_result
        sample_clean = sample_result
        
        # Ensure sequences are the same length for comparison
        min_length = min(len(ref_clean), len(sample_clean))
        if min_length == 0:
            return {
                'success': False,
                'error': "One or both sequences are empty"
            }
        
        # Truncate to same length
        ref_clean = ref_clean[:min_length]
        sample_clean = sample_clean[:min_length]
        
        # Find SNPs
        snps = []
        transitions = 0
        transversions = 0
        
        # Transition pairs (purine to purine, pyrimidine to pyrimidine)
        transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        
        for i in range(min_length):
            ref_base = ref_clean[i]
            sample_base = sample_clean[i]
            
            if ref_base != sample_base and ref_base != 'N' and sample_base != 'N':
                snp_type = 'transition' if (ref_base, sample_base) in transition_pairs else 'transversion'
                
                snps.append({
                    'position': i + 1,
                    'reference': ref_base,
                    'sample': sample_base,
                    'type': snp_type
                })
                
                if snp_type == 'transition':
                    transitions += 1
                else:
                    transversions += 1
        
        # Calculate similarity
        matches = min_length - len(snps)
        similarity = (matches / min_length) * 100 if min_length > 0 else 0
        
        # Calculate Ti/Tv ratio
        ti_tv_ratio = transitions / transversions if transversions > 0 else float('inf') if transitions > 0 else 0
        
        # SNP type breakdown
        snp_types = {
            'A→G': 0, 'G→A': 0, 'C→T': 0, 'T→C': 0,  # Transitions
            'A→T': 0, 'T→A': 0, 'A→C': 0, 'C→A': 0,  # Transversions
            'G→T': 0, 'T→G': 0, 'G→C': 0, 'C→G': 0   # Transversions
        }
        
        for snp in snps:
            mutation = f"{snp['reference']}→{snp['sample']}"
            if mutation in snp_types:
                snp_types[mutation] += 1
        
        return {
            'success': True,
            'data': {
                'total_snps': len(snps),
                'transitions': transitions,
                'transversions': transversions,
                'ti_tv_ratio': round(ti_tv_ratio, 2) if ti_tv_ratio != float('inf') else 0,
                'similarity': round(similarity, 1),
                'sequence_length': min_length,
                'snps': snps,
                'snp_types': snp_types,
                'reference_length': len(ref_clean),
                'sample_length': len(sample_clean)
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f"SNP analysis error: {str(e)}"
        }

def parse_fasta(content):
    """Parse FASTA format content"""
    lines = content.strip().split('\n')
    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    return sequence

def parse_genbank(content):
    """Parse GenBank format content"""
    lines = content.strip().split('\n')
    sequence = ''
    in_sequence = False
    
    for line in lines:
        if line.startswith('ORIGIN'):
            in_sequence = True
            continue
        if line.startswith('//'):
            break
        if in_sequence:
            # Remove line numbers and spaces
            clean_line = re.sub(r'^\s*\d+\s*', '', line).replace(' ', '')
            sequence += clean_line
    
    return sequence

def find_stop_codons(sequence):
    """Find stop codons in a DNA sequence"""
    try:
        is_valid, clean_seq = validate_dna_sequence(sequence)
        if not is_valid:
            return {'success': False, 'error': clean_seq}
        
        stop_codons = ['TAA', 'TAG', 'TGA']
        stop_codon_positions = []
        
        # Check all reading frames
        for frame in range(3):
            for i in range(frame, len(clean_seq) - 2, 3):
                codon = clean_seq[i:i+3]
                if codon in stop_codons:
                    stop_codon_positions.append({
                        'position': i + 1,  # 1-based position
                        'codon': codon,
                        'frame': frame + 1
                    })
        
        return {
            'success': True,
            'data': {
                'stop_codons': len(stop_codon_positions),
                'positions': stop_codon_positions
            }
        }
        
    except Exception as e:
        return {'success': False, 'error': f"Stop codon analysis error: {str(e)}"}

def analyze_stop_codons(reference_seq, sample_seq):
    """Analyze stop codons in both reference and sample sequences"""
    try:
        ref_result = find_stop_codons(reference_seq)
        sample_result = find_stop_codons(sample_seq)
        
        if not ref_result['success']:
            return ref_result
        if not sample_result['success']:
            return sample_result
        
        return {
            'success': True,
            'data': {
                'reference_stop_codons': ref_result['data']['stop_codons'],
                'sample_stop_codons': sample_result['data']['stop_codons'],
                'reference_positions': ref_result['data']['positions'],
                'sample_positions': sample_result['data']['positions']
            }
        }
        
    except Exception as e:
        return {'success': False, 'error': f"Stop codon comparison error: {str(e)}"}

def analyze_coding_impact(reference_seq, sample_seq):
    """Analyze coding impact of SNPs (synonymous/non-synonymous)"""
    try:
        # First get the SNPs
        snp_result = find_snps(reference_seq, sample_seq)
        if not snp_result['success']:
            return snp_result
        
        snps = snp_result['data']['snps']
        
        # Genetic code table
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
        
        is_valid_ref, ref_clean = validate_dna_sequence(reference_seq)
        is_valid_sample, sample_clean = validate_dna_sequence(sample_seq)
        
        if not is_valid_ref or not is_valid_sample:
            return {'success': False, 'error': 'Invalid sequences'}
        
        synonymous = 0
        non_synonymous = 0
        stop_gained = 0
        start_lost = 0
        
        # Detailed mutation information
        mutations = []
        synonymous_mutations = []
        non_synonymous_mutations = []
        stop_gained_mutations = []
        start_lost_mutations = []
        
        for snp in snps:
            pos = snp['position'] - 1  # Convert to 0-based
            
            # Find which codon this position belongs to
            codon_start = (pos // 3) * 3
            codon_pos = pos % 3
            
            if codon_start + 2 < len(ref_clean) and codon_start + 2 < len(sample_clean):
                # Get reference and sample codons
                ref_codon = ref_clean[codon_start:codon_start+3]
                sample_codon = sample_clean[codon_start:codon_start+3]
                
                # Get amino acids
                ref_aa = genetic_code.get(ref_codon, 'X')
                sample_aa = genetic_code.get(sample_codon, 'X')
                
                # Create mutation info
                mutation_info = {
                    'position': pos + 1,  # 1-based position
                    'codon_start': codon_start + 1,  # 1-based codon start
                    'codon_position': codon_pos + 1,  # Position within codon (1-3)
                    'ref_codon': ref_codon,
                    'sample_codon': sample_codon,
                    'ref_aa': ref_aa,
                    'sample_aa': sample_aa,
                    'ref_nucleotide': snp['reference'],
                    'sample_nucleotide': snp['sample'],
                    'change': f"{ref_codon}→{sample_codon}"
                }
                
                # Classify the mutation
                if ref_aa == sample_aa:
                    synonymous += 1
                    mutation_info['type'] = 'synonymous'
                    synonymous_mutations.append(mutation_info)
                else:
                    non_synonymous += 1
                    mutation_info['type'] = 'non_synonymous'
                    non_synonymous_mutations.append(mutation_info)
                    
                    # Check for stop gained
                    if ref_aa != '*' and sample_aa == '*':
                        stop_gained += 1
                        mutation_info['special_type'] = 'stop_gained'
                        stop_gained_mutations.append(mutation_info)
                    
                    # Check for start lost
                    if ref_codon == 'ATG' and sample_codon != 'ATG':
                        start_lost += 1
                        mutation_info['special_type'] = 'start_lost'
                        start_lost_mutations.append(mutation_info)
                
                mutations.append(mutation_info)
        
        return {
            'success': True,
            'data': {
                'synonymous': synonymous,
                'non_synonymous': non_synonymous,
                'stop_gained': stop_gained,
                'start_lost': start_lost,
                'total_analyzed': len(snps),
                'mutations': mutations,
                'synonymous_mutations': synonymous_mutations,
                'non_synonymous_mutations': non_synonymous_mutations,
                'stop_gained_mutations': stop_gained_mutations,
                'start_lost_mutations': start_lost_mutations
            }
        }
        
    except Exception as e:
        return {'success': False, 'error': f"Coding impact analysis error: {str(e)}"}

def main():
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
            seq_type = sys.argv[3] if len(sys.argv) > 3 else 'unknown'
            
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
            result = analyze_sequence(sequence, seq_type)
            print(json.dumps(result))
            
        elif operation == 'analyze':
            # Analyze SNPs between two sequences
            if len(sys.argv) < 4:
                print(json.dumps({
                    'success': False,
                    'error': 'Insufficient arguments for SNP analysis'
                }))
                return
            
            reference_seq = sys.argv[2]
            sample_seq = sys.argv[3]
            
            # Parse sequences if they're in file formats
            if reference_seq.startswith('>'):
                reference_seq = parse_fasta(reference_seq)
            elif 'ORIGIN' in reference_seq:
                reference_seq = parse_genbank(reference_seq)
            
            if sample_seq.startswith('>'):
                sample_seq = parse_fasta(sample_seq)
            elif 'ORIGIN' in sample_seq:
                sample_seq = parse_genbank(sample_seq)
            
            # Find SNPs
            result = find_snps(reference_seq, sample_seq)
            print(json.dumps(result))
            
        elif operation == 'stop_codons':
            # Analyze stop codons in both sequences
            if len(sys.argv) < 4:
                print(json.dumps({
                    'success': False,
                    'error': 'Insufficient arguments for stop codon analysis'
                }))
                return
            
            reference_seq = sys.argv[2]
            sample_seq = sys.argv[3]
            
            # Parse sequences if they're in file formats
            if reference_seq.startswith('>'):
                reference_seq = parse_fasta(reference_seq)
            elif 'ORIGIN' in reference_seq:
                reference_seq = parse_genbank(reference_seq)
            
            if sample_seq.startswith('>'):
                sample_seq = parse_fasta(sample_seq)
            elif 'ORIGIN' in sample_seq:
                sample_seq = parse_genbank(sample_seq)
            
            # Analyze stop codons
            result = analyze_stop_codons(reference_seq, sample_seq)
            print(json.dumps(result))
            
        elif operation == 'coding_impact':
            # Analyze coding impact (synonymous/non-synonymous)
            if len(sys.argv) < 4:
                print(json.dumps({
                    'success': False,
                    'error': 'Insufficient arguments for coding impact analysis'
                }))
                return
            
            reference_seq = sys.argv[2]
            sample_seq = sys.argv[3]
            
            # Parse sequences if they're in file formats
            if reference_seq.startswith('>'):
                reference_seq = parse_fasta(reference_seq)
            elif 'ORIGIN' in reference_seq:
                reference_seq = parse_genbank(reference_seq)
            
            if sample_seq.startswith('>'):
                sample_seq = parse_fasta(sample_seq)
            elif 'ORIGIN' in sample_seq:
                sample_seq = parse_genbank(sample_seq)
            
            # Analyze coding impact
            result = analyze_coding_impact(reference_seq, sample_seq)
            print(json.dumps(result))
            
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