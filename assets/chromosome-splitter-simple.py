#!/usr/bin/env python3
"""

"""

import argparse
import json
import sys
import os
import re
from typing import Dict, List

# Try to import BioPython for better parsing
try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

def parse_fasta_simple(fasta_input) -> Dict[str, str]:
    """Simple FASTA parsing - fast and memory efficient."""
    sequences = {}
    
    if isinstance(fasta_input, str) and os.path.exists(fasta_input):
        # File path - most memory efficient
        if BIOPYTHON_AVAILABLE:
            print("Using BioPython file parsing", file=sys.stderr)
            for record in SeqIO.parse(fasta_input, "fasta"):
                # Just like script: record.id.split()[0]
                chrom_name = record.id.split()[0]
                sequences[chrom_name] = str(record.seq)
        else:
            print("Using manual file parsing", file=sys.stderr)
            with open(fasta_input, 'r') as f:
                current_header = None
                current_sequence = []
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_header and current_sequence:
                            sequences[current_header] = ''.join(current_sequence)
                        
                        # Start new sequence - just like script
                        current_header = line[1:].split()[0]  # Remove '>' and split
                        current_sequence = []
                    else:
                        current_sequence.append(line.upper())
                
                # Save last sequence
                if current_header and current_sequence:
                    sequences[current_header] = ''.join(current_sequence)
    else:
        # Content string
        if BIOPYTHON_AVAILABLE:
            print("Using BioPython string parsing", file=sys.stderr)
            from io import StringIO
            fasta_io = StringIO(fasta_input)
            for record in SeqIO.parse(fasta_io, "fasta"):
                chrom_name = record.id.split()[0]
                sequences[chrom_name] = str(record.seq)
        else:
            print("Using manual string parsing", file=sys.stderr)
            lines = fasta_input.strip().split('\n')
            current_header = None
            current_sequence = []
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    if current_header and current_sequence:
                        sequences[current_header] = ''.join(current_sequence)
                    
                    current_header = line[1:].split()[0]
                    current_sequence = []
                else:
                    current_sequence.append(line.upper())
            
            if current_header and current_sequence:
                sequences[current_header] = ''.join(current_sequence)
    
    return sequences

def classify_chromosome_type(chr_name: str) -> str:
    """Simple chromosome classification."""
    chr_name_lower = chr_name.lower()
    
    if 'chr' in chr_name_lower:
        if any(x in chr_name_lower for x in ['x', 'y', 'z', 'w']):
            return 'sex_chromosome'
        elif any(x in chr_name_lower for x in ['mt', 'mito', 'chloro', 'cp']):
            return 'organellar'
        else:
            return 'autosome'
    elif 'scaffold' in chr_name_lower:
        return 'scaffold'
    elif 'contig' in chr_name_lower:
        return 'contig'
    elif 'unplaced' in chr_name_lower or 'random' in chr_name_lower:
        return 'unplaced'
    else:
        return 'chromosome'

def split_chromosomes_simple(fasta_input) -> Dict:
    """Simple chromosome splitting - just like the working script!"""
    
    try:
        print("Starting simple chromosome splitting", file=sys.stderr)
        
        # Parse FASTA - simple and fast
        sequences = parse_fasta_simple(fasta_input)
        if not sequences:
            return {
                'success': False,
                'error': 'No sequences found in FASTA file'
            }
        
        print(f"Found {len(sequences)} sequences", file=sys.stderr)
        
        # Process each sequence as a chromosome (just like script!)
        chromosome_list = []
        chromosome_files = {}
        
        for chr_name, sequence in sequences.items():
            # Classify chromosome type
            chr_type = classify_chromosome_type(chr_name)
            
            # Generate filename (clean name)
            clean_name = re.sub(r'[^\w\-_.]', '_', chr_name)
            filename = f"{clean_name}.fa"
            
            # Create FASTA content (just like script!)
            fasta_content = f">{chr_name}\n"
            # Add line breaks every 80 characters
            for i in range(0, len(sequence), 80):
                fasta_content += sequence[i:i+80] + "\n"
            
            chromosome_files[filename] = fasta_content
            
            chromosome_list.append({
                'name': chr_name,
                'filename': filename,
                'length': len(sequence),
                'type': chr_type,
                'size_mb': round(len(sequence) / 1_000_000, 2),
                'coordinates': f"1-{len(sequence)}",
                'original_header': chr_name
            })
            
            print(f"✅ Processed: {chr_name} ({len(sequence):,} bp) - {chr_type}", file=sys.stderr)
        
        # Sort chromosomes for better display
        type_order = {'autosome': 1, 'sex_chromosome': 2, 'organellar': 3, 'chromosome': 4, 'scaffold': 5, 'contig': 6, 'unplaced': 7}
        chromosome_list.sort(key=lambda x: (type_order.get(x['type'], 8), x['name']))
        
        print(f"🎉 Successfully split {len(sequences)} chromosomes!", file=sys.stderr)
        
        return {
            'success': True,
            'chromosomes': chromosome_list,
            'files': chromosome_files,
            'statistics': {
                'total_sequences': len(sequences),
                'total_bases': sum(len(seq) for seq in sequences.values()),
                'chromosomes_found': len(sequences),
                'files_created': len(chromosome_files)
            },
            'summary': {
                'total_chromosomes': len(sequences),
                'total_size_mb': round(sum(len(seq) for seq in sequences.values()) / 1_000_000, 2),
                'file_count': len(chromosome_files)
            }
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': f'Processing failed: {str(e)}'
        }

def main():
    parser = argparse.ArgumentParser(description='Simple Chromosome Splitter - Fast FASTA splitting')
    
    # FASTA input options
    parser.add_argument('--fasta-content', help='FASTA file content (direct)')
    parser.add_argument('--fasta-file', help='Path to FASTA file')
    parser.add_argument('--output-dir', help='Output directory for chromosome files')
    
    args = parser.parse_args()
    
    # Get FASTA input
    fasta_input = None
    if args.fasta_file and os.path.exists(args.fasta_file):
        fasta_input = args.fasta_file
        print(f"Using FASTA file: {args.fasta_file}", file=sys.stderr)
    elif args.fasta_content:
        fasta_input = args.fasta_content
        print("Using FASTA content string", file=sys.stderr)
    elif args.fasta_file:
        print(json.dumps({
            'success': False,
            'error': f'FASTA file not found: {args.fasta_file}'
        }))
        return
    else:
        print(json.dumps({
            'success': False,
            'error': 'Either --fasta-content or --fasta-file must be provided'
        }))
        return
    
    # Split chromosomes
    result = split_chromosomes_simple(fasta_input)
    
    # Save files if output directory specified
    if args.output_dir and result['success']:
        try:
            os.makedirs(args.output_dir, exist_ok=True)
            
            for filename, content in result['files'].items():
                filepath = os.path.join(args.output_dir, filename)
                with open(filepath, 'w') as f:
                    f.write(content)
                print(f"Saved: {filepath}", file=sys.stderr)
                
        except Exception as e:
            result['error'] = f'Failed to save files: {str(e)}'
            result['success'] = False
    
    # Output result as JSON
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()