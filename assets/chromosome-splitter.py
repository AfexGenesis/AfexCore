#!/usr/bin/env python3
"""
Chromosome Splitter - Advanced genome splitting tool
Uses FASTA genome files and GFF3 annotations to split chromosomes accurately
"""

import argparse
import json
import sys
import re
import os
import tempfile
import zipfile
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

class ChromosomeSplitter:
    def __init__(self):
        """Initialize the Chromosome Splitter."""
        self.chromosomes = {}
        self.sequences = {}
        self.statistics = {
            'total_sequences': 0,
            'total_bases': 0,
            'chromosomes_found': 0,
            'scaffolds_found': 0,
            'files_created': 0
        }
    
    def parse_fasta(self, fasta_content: str) -> Dict[str, str]:
        """Parse FASTA content and return sequences."""
        sequences = {}
        current_header = None
        current_sequence = []
        
        lines = fasta_content.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence
                if current_header and current_sequence:
                    sequences[current_header] = ''.join(current_sequence)
                
                # Start new sequence
                current_header = line[1:].strip()  # Remove '>'
                current_sequence = []
            else:
                # Add to current sequence
                current_sequence.append(line.upper())
        
        # Save last sequence
        if current_header and current_sequence:
            sequences[current_header] = ''.join(current_sequence)
        
        return sequences
    
    def parse_gff3(self, gff_content: str) -> List[Dict]:
        """Parse GFF3 content and extract chromosome features."""
        features = []
        
        lines = gff_content.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            # Split GFF3 fields
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip()
            
            feature = {
                'seqid': seqid,
                'source': source,
                'type': feature_type.lower(),
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attr_dict
            }
            
            features.append(feature)
        
        return features
    
    def extract_chromosomes_from_gff(self, gff_features: List[Dict]) -> Dict[str, Dict]:
        """Extract chromosome information from GFF3 features."""
        chromosomes = {}
        
        for feature in gff_features:
            feature_type = feature['type']
            seqid = feature['seqid']
            
            # Look for chromosome, scaffold, or contig features
            if feature_type in ['chromosome', 'scaffold', 'contig', 'supercontig']:
                
                # Get chromosome name from attributes
                chr_name = None
                attributes = feature['attributes']
                
                # Try different attribute names
                for attr_key in ['Name', 'ID', 'chromosome', 'chr']:
                    if attr_key in attributes:
                        chr_name = attributes[attr_key]
                        break
                
                # If no name found, use seqid
                if not chr_name:
                    chr_name = seqid
                
                chromosomes[chr_name] = {
                    'seqid': seqid,
                    'type': feature_type,
                    'start': feature['start'],
                    'end': feature['end'],
                    'length': feature['end'] - feature['start'] + 1,
                    'attributes': attributes
                }
        
        return chromosomes
    
    def map_chromosomes_to_sequences(self, chromosomes: Dict, sequences: Dict) -> Dict:
        """Map chromosome information to FASTA sequences."""
        mapped_chromosomes = {}
        
        for chr_name, chr_info in chromosomes.items():
            seqid = chr_info['seqid']
            
            # Try to find matching sequence
            sequence = None
            matched_header = None
            
            # Direct match
            if seqid in sequences:
                sequence = sequences[seqid]
                matched_header = seqid
            else:
                # Try partial matches
                for header, seq in sequences.items():
                    if seqid in header or header in seqid:
                        sequence = seq
                        matched_header = header
                        break
                
                # Try case-insensitive match
                if not sequence:
                    for header, seq in sequences.items():
                        if seqid.lower() in header.lower() or header.lower() in seqid.lower():
                            sequence = seq
                            matched_header = header
                            break
            
            if sequence:
                # Extract chromosome sequence based on coordinates
                start = max(0, chr_info['start'] - 1)  # Convert to 0-based
                end = min(len(sequence), chr_info['end'])
                
                chr_sequence = sequence[start:end] if start < end else sequence
                
                mapped_chromosomes[chr_name] = {
                    'sequence': chr_sequence,
                    'length': len(chr_sequence),
                    'type': chr_info['type'],
                    'original_header': matched_header,
                    'coordinates': f"{chr_info['start']}-{chr_info['end']}",
                    'attributes': chr_info['attributes']
                }
        
        return mapped_chromosomes
    
    def classify_chromosome_type(self, chr_name: str, chr_info: Dict) -> str:
        """Classify chromosome type for better organization."""
        chr_name_lower = chr_name.lower()
        chr_type = chr_info.get('type', '').lower()
        
        # Check for chromosome patterns
        if 'chr' in chr_name_lower or chr_type == 'chromosome':
            if any(x in chr_name_lower for x in ['x', 'y', 'z', 'w']):
                return 'sex_chromosome'
            elif any(x in chr_name_lower for x in ['mt', 'mito', 'chloro', 'cp']):
                return 'organellar'
            else:
                return 'autosome'
        
        # Check for scaffold/contig patterns
        elif 'scaffold' in chr_name_lower or chr_type == 'scaffold':
            return 'scaffold'
        elif 'contig' in chr_name_lower or chr_type == 'contig':
            return 'contig'
        else:
            return 'unknown'
    
    def generate_output_filename(self, chr_name: str, chr_info: Dict) -> str:
        """Generate appropriate filename for chromosome."""
        chr_type = self.classify_chromosome_type(chr_name, chr_info)
        
        # Clean chromosome name for filename
        clean_name = re.sub(r'[^\w\-_.]', '_', chr_name)
        
        # Add appropriate prefix
        if chr_type == 'autosome':
            return f"Chr_{clean_name}.fasta"
        elif chr_type == 'sex_chromosome':
            return f"Chr_{clean_name}_sex.fasta"
        elif chr_type == 'organellar':
            return f"Organelle_{clean_name}.fasta"
        elif chr_type == 'scaffold':
            return f"Scaffold_{clean_name}.fasta"
        elif chr_type == 'contig':
            return f"Contig_{clean_name}.fasta"
        else:
            return f"{clean_name}.fasta"
    
    def create_fasta_content(self, chr_name: str, chr_info: Dict) -> str:
        """Create FASTA content for a chromosome."""
        sequence = chr_info['sequence']
        length = chr_info['length']
        chr_type = chr_info['type']
        coordinates = chr_info.get('coordinates', '')
        
        # Create descriptive header
        header = f">{chr_name}"
        if coordinates:
            header += f" coordinates:{coordinates}"
        header += f" length:{length} type:{chr_type}"
        
        # Format sequence with line breaks (80 characters per line)
        formatted_sequence = ""
        for i in range(0, len(sequence), 80):
            formatted_sequence += sequence[i:i+80] + "\n"
        
        return header + "\n" + formatted_sequence
    
    def calculate_statistics(self, sequences: Dict, chromosomes: Dict):
        """Calculate processing statistics."""
        self.statistics['total_sequences'] = len(sequences)
        self.statistics['total_bases'] = sum(len(seq) for seq in sequences.values())
        self.statistics['chromosomes_found'] = len(chromosomes)
        
        # Count different types
        scaffold_count = 0
        for chr_info in chromosomes.values():
            chr_type = self.classify_chromosome_type('', chr_info)
            if chr_type in ['scaffold', 'contig']:
                scaffold_count += 1
        
        self.statistics['scaffolds_found'] = scaffold_count
        self.statistics['files_created'] = len(chromosomes)
    
    def split_chromosomes(self, fasta_content: str, gff_content: str) -> Dict:
        """Main function to split chromosomes using FASTA and GFF3."""
        
        try:
            # Parse input files
            sequences = self.parse_fasta(fasta_content)
            if not sequences:
                return {
                    'success': False,
                    'error': 'No sequences found in FASTA file'
                }
            
            gff_features = self.parse_gff3(gff_content)
            if not gff_features:
                return {
                    'success': False,
                    'error': 'No features found in GFF3 file'
                }
            
            # Extract chromosome information from GFF3
            chromosomes = self.extract_chromosomes_from_gff(gff_features)
            if not chromosomes:
                return {
                    'success': False,
                    'error': 'No chromosome features found in GFF3 file'
                }
            
            # Map chromosomes to sequences
            mapped_chromosomes = self.map_chromosomes_to_sequences(chromosomes, sequences)
            if not mapped_chromosomes:
                return {
                    'success': False,
                    'error': 'Could not map GFF3 chromosomes to FASTA sequences'
                }
            
            # Calculate statistics
            self.calculate_statistics(sequences, mapped_chromosomes)
            
            # Prepare results
            chromosome_files = {}
            chromosome_list = []
            
            for chr_name, chr_info in mapped_chromosomes.items():
                filename = self.generate_output_filename(chr_name, chr_info)
                fasta_content = self.create_fasta_content(chr_name, chr_info)
                
                chromosome_files[filename] = fasta_content
                
                chromosome_list.append({
                    'name': chr_name,
                    'filename': filename,
                    'length': chr_info['length'],
                    'type': self.classify_chromosome_type(chr_name, chr_info),
                    'size_mb': round(chr_info['length'] / 1_000_000, 2),
                    'coordinates': chr_info.get('coordinates', ''),
                    'original_header': chr_info.get('original_header', '')
                })
            
            # Sort chromosome list by type and name
            chromosome_list.sort(key=lambda x: (x['type'], x['name']))
            
            return {
                'success': True,
                'chromosomes': chromosome_list,
                'files': chromosome_files,
                'statistics': self.statistics,
                'summary': {
                    'total_chromosomes': len(mapped_chromosomes),
                    'total_size_mb': round(sum(info['length'] for info in mapped_chromosomes.values()) / 1_000_000, 2),
                    'file_count': len(chromosome_files)
                }
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f'Processing failed: {str(e)}'
            }

def main():
    parser = argparse.ArgumentParser(description='Chromosome Splitter - Split genome using FASTA and GFF3')
    parser.add_argument('--fasta-content', required=True, help='FASTA file content')
    parser.add_argument('--gff-content', required=True, help='GFF3 file content')
    parser.add_argument('--output-dir', help='Output directory for chromosome files')
    
    args = parser.parse_args()
    
    splitter = ChromosomeSplitter()
    
    # Split chromosomes
    result = splitter.split_chromosomes(args.fasta_content, args.gff_content)
    
    # Save files if output directory specified
    if result['success'] and args.output_dir and 'files' in result:
        try:
            os.makedirs(args.output_dir, exist_ok=True)
            
            for filename, content in result['files'].items():
                filepath = os.path.join(args.output_dir, filename)
                with open(filepath, 'w') as f:
                    f.write(content)
            
            result['output_directory'] = args.output_dir
            result['files_saved'] = len(result['files'])
            
        except Exception as e:
            result['warning'] = f'Files created but could not save to disk: {str(e)}'
    
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()