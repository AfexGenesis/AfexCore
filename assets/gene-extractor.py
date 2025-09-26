#!/usr/bin/env python3
"""
AFEX Gene Extractor - FAST Two-Stage Approach
Stage 1: Quick scan for gene info (FAST)
Stage 2: Extract sequences on-demand (INSTANT)
Based on afexgenesnatcher.py for maximum speed
"""

import os
import sys
import json
import time
import argparse
from io import StringIO

try:
    from BCBio import GFF
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError as e:
    print(json.dumps({
        'success': False,
        'error': f'Missing required library: {str(e)}. Please install: pip install biopython bcbio-gff'
    }))
    sys.exit(1)


class AFEXFastGeneExtractor:
    """FAST gene extractor using afexgenesnatcher.py approach"""
    
    def __init__(self):
        self.genome_dict = None
        self.genome_file_path = None
    
    def quick_gene_scan(self, genome_data, gff_data, file_name='genome', page=1, genes_per_page=500, load_all=False):
        """STAGE 1: Quick scan to get gene info WITHOUT extracting sequences"""
        try:
            start_time = time.time()
            sys.stderr.write("🚀 FAST GENE SCANNER - Stage 1: Quick Info Scan\n")
            sys.stderr.write("=" * 50 + "\n")
            sys.stderr.flush()
            
            # Store genome reference for later extraction
            if os.path.exists(genome_data):
                self.genome_file_path = genome_data
                sys.stderr.write(f"📁 Genome file: {genome_data}\n")
            else:
                sys.stderr.write("📝 Genome from content data\n")
            sys.stderr.flush()
            
            # Parse GFF ONLY (no genome loading yet - FAST!)
            sys.stderr.write("📋 Quick GFF scan (no sequence extraction)...\n")
            sys.stderr.flush()
            
            gene_info_list = []
            skipped = 0
            processed_features = 0
            
            # Create temp GFF file if needed
            if not os.path.exists(gff_data):
                import tempfile
                gff_temp = tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False, encoding='utf-8')
                gff_temp.write(gff_data)
                gff_temp.close()
                gff_file_path = gff_temp.name
            else:
                gff_file_path = gff_data
            
            # FAST GFF parsing - INFO ONLY (like afexgenesnatcher.py but no sequence extraction)
            with open(gff_file_path, 'r') as gff_handle:
                for rec in GFF.parse(gff_handle):  # NO base_dict = FAST!
                    sys.stderr.write(f"🧬 Scanning: {rec.id}\n")
                    sys.stderr.flush()
                    
                    for feature in rec.features:
                        processed_features += 1
                        
                        if feature.type in ["gene", "mRNA", "CDS"]:
                            try:
                                # Extract gene INFO only (no sequences)
                                qualifiers = feature.qualifiers
                                
                                # Debug: Show available attributes for first few genes
                                if len(gene_info_list) < 3:
                                    sys.stderr.write(f"🔍 Gene {len(gene_info_list)+1} attributes: {list(qualifiers.keys())}\n")
                                    for attr in ["ID", "Name", "gene", "gene_id", "locus_tag"]:
                                        if attr in qualifiers:
                                            sys.stderr.write(f"   {attr}: {qualifiers[attr]}\n")
                                    sys.stderr.flush()
                                
                                # Try multiple common GFF attributes for gene ID
                                gene_id = None
                                for id_attr in ["ID", "gene_id", "locus_tag", "Name"]:
                                    if id_attr in qualifiers and qualifiers[id_attr]:
                                        gene_id = qualifiers[id_attr][0].replace(":", "_").replace(";", "_")
                                        break
                                if not gene_id:
                                    gene_id = f"unknown_gene_{processed_features}"
                                
                                # Try multiple common GFF attributes for gene name
                                gene_name = None
                                for name_attr in ["Name", "gene", "gene_name", "locus_tag", "ID"]:
                                    if name_attr in qualifiers and qualifiers[name_attr]:
                                        gene_name = qualifiers[name_attr][0]
                                        break
                                if not gene_name:
                                    gene_name = gene_id
                                
                                # Try multiple common GFF attributes for product/function
                                product = None
                                for prod_attr in ["product", "function", "description", "note", "annotation"]:
                                    if prod_attr in qualifiers and qualifiers[prod_attr]:
                                        product = qualifiers[prod_attr][0]
                                        break
                                if not product:
                                    product = "Unknown protein"
                                
                                # Debug: Show what was selected for first few genes
                                if len(gene_info_list) < 3:
                                    sys.stderr.write(f"   Selected ID: {gene_id}\n")
                                    sys.stderr.write(f"   Selected Name: {gene_name}\n")
                                    sys.stderr.write(f"   Selected Product: {product}\n")
                                    sys.stderr.flush()
                                
                                start = int(feature.location.start) + 1  # 1-based
                                end = int(feature.location.end)
                                strand = "+" if feature.location.strand == 1 else "-"
                                length = end - start + 1
                                
                                gene_info = {
                                    'id': gene_id,
                                    'name': gene_name,
                                    'type': feature.type,
                                    'chromosome': rec.id,
                                    'start': start,
                                    'end': end,
                                    'strand': strand,
                                    'length': length,
                                    'product': product,
                                    'description': f"{gene_name} | {rec.id}:{start}-{end}({strand}) | {product}",
                                    'sequence': None  # Will be extracted on-demand
                                }
                                
                                gene_info_list.append(gene_info)
                                
                            except Exception as e:
                                skipped += 1
                                if skipped <= 5:  # Show first 5 errors only
                                    sys.stderr.write(f"⚠️ Skipped: {str(e)}\n")
                                    sys.stderr.flush()
            
            # Clean up temp file
            if not os.path.exists(gff_data):
                try:
                    os.remove(gff_file_path)
                except:
                    pass
            
            elapsed = round(time.time() - start_time, 2)
            
            sys.stderr.write("=" * 50 + "\n")
            sys.stderr.write(f"⚡ QUICK SCAN COMPLETE!\n")
            sys.stderr.write(f"✅ Found {len(gene_info_list):,} genes in {elapsed}s\n")
            sys.stderr.write(f"❌ Skipped: {skipped:,}\n")
            sys.stderr.write("🎯 Ready for on-demand sequence extraction!\n")
            sys.stderr.write("=" * 50 + "\n")
            sys.stderr.flush()
            
            if not gene_info_list:
                return {
                    'success': False,
                    'error': f'No genes found. Processed {processed_features} features.'
                }
            
            # Calculate basic statistics (fast)
            gene_types = {}
            chr_counts = {}
            lengths = []
            
            for gene in gene_info_list:
                gene_types[gene['type']] = gene_types.get(gene['type'], 0) + 1
                chr_counts[gene['chromosome']] = chr_counts.get(gene['chromosome'], 0) + 1
                lengths.append(gene['length'])
            
            statistics = {
                'total_genes': len(gene_info_list),
                'average_length': round(sum(lengths) / len(lengths), 2),
                'longest_gene': max(lengths),
                'shortest_gene': min(lengths),
                'gene_types': gene_types,
                'chromosome_distribution': chr_counts,
                'scan_time': elapsed
            }
            
            # Handle pagination vs load all
            total_genes = len(gene_info_list)
            
            if load_all:
                # Load ALL genes at once - let's see if it crashes!
                display_genes = gene_info_list
                total_pages = 1
                start_idx = 0
                end_idx = total_genes
                sys.stderr.write(f"🔥 Loading ALL {total_genes:,} genes at once!\n")
                sys.stderr.flush()
            else:
                # Normal pagination
                total_pages = (total_genes + genes_per_page - 1) // genes_per_page
                start_idx = (page - 1) * genes_per_page
                end_idx = min(start_idx + genes_per_page, total_genes)
                display_genes = gene_info_list[start_idx:end_idx]
                sys.stderr.write(f"📄 Page {page}/{total_pages}: Showing genes {start_idx+1}-{end_idx} of {total_genes}\n")
                sys.stderr.flush()
            
            return {
                'success': True,
                'genes': display_genes,
                'statistics': statistics,
                'method': 'FAST_GFF_scan',
                'total_genes_found': len(gene_info_list),
                'genes_skipped': skipped,
                'processing_time': elapsed,
                'pagination': {
                    'total_genes': total_genes,
                    'genes_per_page': genes_per_page,
                    'total_pages': total_pages,
                    'current_page': page,
                    'showing_genes': len(display_genes),
                    'start_gene': start_idx + 1,
                    'end_gene': end_idx,
                    'has_next': page < total_pages,
                    'has_prev': page > 1
                },
                'note': f'FAST SCAN - Found {len(gene_info_list):,} genes in {elapsed}s. Sequences available on-demand.',
                'extraction_ready': True  # Flag for on-demand extraction
            }
            
        except Exception as e:
            sys.stderr.write(f"❌ Error during quick scan: {str(e)}\n")
            sys.stderr.flush()
            return {
                'success': False,
                'error': f'Quick scan failed: {str(e)}'
            }
    
    def extract_single_gene(self, gene_info, genome_data):
        """STAGE 2: Extract single gene sequence on-demand (INSTANT)"""
        try:
            sys.stderr.write(f"🎯 Extracting: {gene_info['id']}\n")
            sys.stderr.flush()
            
            # Load genome if not already loaded
            if self.genome_dict is None:
                sys.stderr.write("📖 Loading genome for extraction...\n")
                sys.stderr.flush()
                
                if self.genome_file_path and os.path.exists(self.genome_file_path):
                    # Use afexgenesnatcher.py approach - direct load
                    self.genome_dict = SeqIO.to_dict(SeqIO.parse(self.genome_file_path, "fasta"))
                else:
                    # Load from content
                    genome_io = StringIO(genome_data)
                    self.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_io, "fasta"))
                
                sys.stderr.write(f"✅ Loaded {len(self.genome_dict)} sequences\n")
                sys.stderr.flush()
            
            # Extract sequence (like afexgenesnatcher.py)
            chromosome = gene_info['chromosome']
            if chromosome not in self.genome_dict:
                return {
                    'success': False,
                    'error': f'Chromosome {chromosome} not found in genome'
                }
            
            chr_seq = self.genome_dict[chromosome].seq
            start = gene_info['start'] - 1  # Convert to 0-based
            end = gene_info['end']
            strand = gene_info['strand']
            
            # Extract sequence
            gene_seq = chr_seq[start:end]
            if strand == "-":
                gene_seq = gene_seq.reverse_complement()
            
            # Calculate GC content
            gc_count = str(gene_seq).upper().count('G') + str(gene_seq).upper().count('C')
            gc_content = round((gc_count / len(gene_seq)) * 100, 2) if len(gene_seq) > 0 else 0
            
            # Return complete gene info with sequence
            complete_gene = gene_info.copy()
            complete_gene['sequence'] = str(gene_seq)
            complete_gene['gc_content'] = gc_content
            
            sys.stderr.write(f"✅ Extracted {gene_info['id']} ({len(gene_seq)} bp)\n")
            sys.stderr.flush()
            
            return {
                'success': True,
                'gene': complete_gene
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f'Failed to extract {gene_info["id"]}: {str(e)}'
            }
    
    def extract_multiple_genes(self, gene_list, genome_data):
        """STAGE 2: Extract multiple genes on-demand"""
        try:
            sys.stderr.write(f"🎯 Bulk extraction: {len(gene_list)} genes\n")
            sys.stderr.flush()
            
            # Load genome if not already loaded
            if self.genome_dict is None:
                sys.stderr.write("📖 Loading genome for bulk extraction...\n")
                sys.stderr.flush()
                
                if self.genome_file_path and os.path.exists(self.genome_file_path):
                    self.genome_dict = SeqIO.to_dict(SeqIO.parse(self.genome_file_path, "fasta"))
                else:
                    genome_io = StringIO(genome_data)
                    self.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_io, "fasta"))
                
                sys.stderr.write(f"✅ Loaded {len(self.genome_dict)} sequences\n")
                sys.stderr.flush()
            
            extracted_genes = []
            failed_genes = []
            
            for i, gene_info in enumerate(gene_list):
                try:
                    # Progress indicator
                    if (i + 1) % 100 == 0:
                        sys.stderr.write(f"   Extracted {i + 1}/{len(gene_list)} genes...\n")
                        sys.stderr.flush()
                    
                    # Extract like afexgenesnatcher.py
                    chromosome = gene_info['chromosome']
                    if chromosome not in self.genome_dict:
                        failed_genes.append(gene_info['id'])
                        continue
                    
                    chr_seq = self.genome_dict[chromosome].seq
                    start = gene_info['start'] - 1
                    end = gene_info['end']
                    strand = gene_info['strand']
                    
                    gene_seq = chr_seq[start:end]
                    if strand == "-":
                        gene_seq = gene_seq.reverse_complement()
                    
                    # Calculate GC content
                    gc_count = str(gene_seq).upper().count('G') + str(gene_seq).upper().count('C')
                    gc_content = round((gc_count / len(gene_seq)) * 100, 2) if len(gene_seq) > 0 else 0
                    
                    complete_gene = gene_info.copy()
                    complete_gene['sequence'] = str(gene_seq)
                    complete_gene['gc_content'] = gc_content
                    
                    extracted_genes.append(complete_gene)
                    
                except Exception as e:
                    failed_genes.append(gene_info['id'])
            
            sys.stderr.write(f"✅ Bulk extraction complete: {len(extracted_genes)} success, {len(failed_genes)} failed\n")
            sys.stderr.flush()
            
            return {
                'success': True,
                'genes': extracted_genes,
                'failed_genes': failed_genes,
                'total_extracted': len(extracted_genes)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f'Bulk extraction failed: {str(e)}'
            }


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(description='AFEX Fast Gene Extractor - Two-Stage Approach')
    parser.add_argument('--sequence-data', help='Genome sequence data (FASTA content)')
    parser.add_argument('--sequence-file', help='Path to genome sequence file')
    parser.add_argument('--gff-data', help='GFF annotation data')
    parser.add_argument('--gff-file', help='Path to GFF file')
    parser.add_argument('--file-name', default='genome', help='Original file name for context')
    parser.add_argument('--mode', choices=['scan', 'extract'], default='scan', help='Operation mode')
    parser.add_argument('--gene-id', help='Gene ID for single extraction')
    parser.add_argument('--page', type=int, default=1, help='Page number for pagination')
    parser.add_argument('--genes-per-page', type=int, default=500, help='Genes per page')
    parser.add_argument('--load-all', action='store_true', help='Load ALL genes at once (test for crashes)')
    
    args = parser.parse_args()
    

    
    extractor = AFEXFastGeneExtractor()
    
    try:
        # Get sequence data
        sequence_data = args.sequence_data
        if args.sequence_file and os.path.exists(args.sequence_file):
            sequence_data = args.sequence_file
        
        # Get GFF data
        gff_data = args.gff_data
        if args.gff_file and os.path.exists(args.gff_file):
            gff_data = args.gff_file
        
        if not sequence_data:
            raise ValueError("No genome sequence data provided")
        
        if not gff_data:
            raise ValueError("No GFF data provided - fast extraction requires GFF annotations")
        
        # Execute based on mode
        if args.mode == 'scan':
            # Stage 1: Quick scan with pagination
            result = extractor.quick_gene_scan(sequence_data, gff_data, args.file_name, args.page, args.genes_per_page, args.load_all)
        else:
            # Stage 2: Extract (would need gene info passed somehow)
            raise ValueError("Extract mode not implemented in CLI - use scan mode")
        
        # Output JSON result
        print(json.dumps(result, ensure_ascii=True))
        
    except Exception as e:
        error_result = {
            'success': False,
            'error': f'Operation failed: {str(e)}'
        }
        print(json.dumps(error_result, ensure_ascii=True))
        sys.exit(1)


if __name__ == '__main__':
    main()