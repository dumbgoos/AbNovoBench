#!/usr/bin/env python3
"""
De Novo Sequencing Results Processing Script

This script processes raw de novo sequencing results from various tools and converts them
to standardized summary formats. This is the FIRST step in the AbNovoBench pipeline and
generates the summary CSV files that are required by subsequent metric analysis scripts.

Supported Tools:
    - CasanovoV1, CasanovoV2, CasanovoV3
    - PepNet, pNovo3 
    - pi-HelixNovo, pi-PrimeNovo
    - AdaNovo, ContraNovo
    - DeepNovo, InstaNovo, InstaNovoPlus
    - PointNovo, PGPointNovo
    - SMSNet

Pipeline Position:
    1. [THIS SCRIPT] Process raw tool outputs → Generate summary CSV files
    2. Assembly scripts (alps_fusion.py, assembly_stitch.py) use summary files
    3. Metric scripts (accuracy_metric.py, robustness_metric.py) use summary files

Output Format:
    Generates standardized CSV files with columns for peptide sequences, scores, 
    and metadata for each tool, which serve as input for downstream analysis.
    
    File Structure:
    - Individual tool summaries: {Tool}_summary.csv (e.g., CasanovoV1_summary.csv)
    - Merged summary file: summary_merged.csv (combines all tools into one file)
    
    The summary_merged.csv file contains:
    - Ground truth sequences from original database search
    - Peptide predictions from each tool ({Tool} Peptide columns)
    - Confidence scores from each tool ({Tool} Score columns)
    - Additional metadata (Spectrum Name, Modified Sequence, etc.)

Usage:
    python denovo_process.py --tool CasanovoV1 --input_path /path/to/raw/results --output_path /path/to/summary
"""

import os
import argparse
import pandas as pd
import numpy as np
import yaml
import json
from typing import Dict, Optional, Any
from pyteomics import mgf
from collections import OrderedDict


class DenovoConfigLoader:
    """Configuration loader for Denovo Process parameters."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration loader.
        
        Args:
            config_path: Path to configuration file. If None, uses default path.
        """
        if config_path is None:
            # Default config path is in the src/config directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            # Go up one level from summary to src, then into config
            src_dir = os.path.dirname(script_dir)
            config_path = os.path.join(src_dir, 'config', 'denovo_process.yaml')
        
        self.config_path = config_path
        self.config = self._load_config()
        
        # Get project root directory for resolving relative paths
        self.project_root = self._get_project_root()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            print(f"Loaded configuration from {self.config_path}")
            return config
        except FileNotFoundError:
            print(f"Warning: Configuration file {self.config_path} not found. Using fallback defaults.")
            return self._get_fallback_config()
        except yaml.YAMLError as e:
            print(f"Error parsing YAML config: {e}. Using fallback defaults.")
            return self._get_fallback_config()
    
    def _get_project_root(self) -> str:
        """Get project root directory."""
        # Start from script directory and go up to find project root
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Go up two levels: summary -> src -> project_root
        project_root = os.path.dirname(os.path.dirname(script_dir))
        return project_root
    
    def _resolve_path(self, relative_path: str) -> str:
        """Resolve relative path to absolute path based on project root."""
        if os.path.isabs(relative_path):
            return relative_path
        return os.path.join(self.project_root, relative_path)

    def _get_fallback_config(self) -> Dict[str, Any]:
        """Get fallback configuration if config file is not available."""
        return {
            'tool_names': [
                'AdaNovo', 'CasanovoV1', 'CasanovoV2', 'ContraNovo', 'DeepNovo',
                'InstaNovo', 'PepNet', 'PGPointNovo', 'pi-HelixNovo', 'pi-PrimeNovo',
                'pNovo3', 'PointNovo', 'SMSNet'
            ],
            'default_names': {
                '50ugmAb1': 'mAb1',
                '100ugmAb2': 'mAb2',
                '200ugmAb3': 'mAb3',
                '20230210-mAb4': 'mAb4',
                '20230707-mAb5': 'mAb5',
                '20231203-mAb6': 'mAb6',
                '20231221-mAb7': 'mAb7',
                '20250415-mAb8': 'mAb8'
            },
            'processing_parameters': {
                'match_mass_tol': 0.1,
                'prefix_mass_tol': 0.5
            },
            'default_paths': {
                'output_base': 'processed_results'
            }
        }
    
    @property
    def tool_names(self) -> set:
        """Get supported tool names."""
        return set(self.config.get('tool_names', []))
    
    @property
    def default_names(self) -> Dict[str, str]:
        """Get default antibody name mapping."""
        return self.config.get('default_names', {})
    
    @property
    def processing_parameters(self) -> Dict[str, Any]:
        """Get processing parameters."""
        return self.config.get('processing_parameters', {})
    
    @property
    def vocabulary(self) -> Dict[str, Any]:
        """Get vocabulary configuration."""
        return self.config.get('vocabulary', {})
    
    @property
    def mass_constants(self) -> Dict[str, Any]:
        """Get mass constants."""
        return self.config.get('mass_constants', {})
    
    @property
    def default_paths(self) -> Dict[str, str]:
        """Get default file paths (resolved to absolute paths)."""
        paths = self.config.get('default_paths', {})
        resolved_paths = {}
        for key, path in paths.items():
            resolved_paths[key] = self._resolve_path(path)
        return resolved_paths
    
    @property
    def base_paths(self) -> Dict[str, str]:
        """Get base directory paths (resolved to absolute paths)."""
        paths = self.config.get('base_paths', {})
        resolved_paths = {}
        for key, path in paths.items():
            resolved_paths[key] = self._resolve_path(path)
        return resolved_paths
    
    @property
    def tool_file_patterns(self) -> Dict[str, Dict[str, Any]]:
        """Get tool-specific file patterns."""
        return self.config.get('tool_file_patterns', {})
    
    @property
    def tool_input_dirs(self) -> Dict[str, str]:
        """Get tool-specific input directory mapping."""
        return self.config.get('tool_input_dirs', {})
    
    @property
    def mgf_file_pattern(self) -> str:
        """Get MGF file pattern."""
        return self.config.get('mgf_file_pattern', 'spectrum_{sample_name}_HCD.mgf')
    
    @property
    def mgf_subdir(self) -> str:
        """Get MGF subdirectory."""
        return self.config.get('mgf_subdir', 'process1')


# Global configuration loader
_denovo_config_loader = None

def get_denovo_config_loader(config_path: Optional[str] = None) -> DenovoConfigLoader:
    """Get global configuration loader instance."""
    global _denovo_config_loader
    if _denovo_config_loader is None or config_path is not None:
        _denovo_config_loader = DenovoConfigLoader(config_path)
    return _denovo_config_loader


class DenovoProcessor:
    """Denovo sequencing result processor for various tools."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the Denovo processor.
        
        Args:
            config_path: Path to configuration file
        """
        # Load configuration
        self.config = get_denovo_config_loader(config_path)
        
        # Set processing parameters from config
        proc_params = self.config.processing_parameters
        self.match_mass_tol = proc_params.get('match_mass_tol', 0.1)
        self.prefix_mass_tol = proc_params.get('prefix_mass_tol', 0.5)
        
        # Load vocabulary
        vocab_config = self.config.vocabulary
        self._setup_vocabulary(vocab_config)
        
        # Load mass constants
        mass_config = self.config.mass_constants
        self._setup_mass_constants(mass_config)
    
    def _setup_vocabulary(self, vocab_config: Dict[str, Any]) -> None:
        """Setup vocabulary from configuration."""
        start_vocab = vocab_config.get('start_vocab', ['_PAD', '_GO', '_EOS'])
        amino_acids = vocab_config.get('amino_acids', [])
        
        self._PAD, self._GO, self._EOS = start_vocab
        self.PAD_ID = vocab_config.get('pad_id', 0)
        self.GO_ID = vocab_config.get('go_id', 1)
        self.EOS_ID = vocab_config.get('eos_id', 2)
        
        self.vocab_reverse = start_vocab + amino_acids
        self.vocab = dict([(x, y) for (y, x) in enumerate(self.vocab_reverse)])
        self.vocab_size = len(self.vocab_reverse)
    
    def _setup_mass_constants(self, mass_config: Dict[str, Any]) -> None:
        """Setup mass constants from configuration."""
        self.mass_H = mass_config.get('mass_H', 1.0078)
        self.mass_H2O = mass_config.get('mass_H2O', 18.0106)
        self.mass_NH3 = mass_config.get('mass_NH3', 17.0265)
        self.mass_N_terminus = mass_config.get('mass_N_terminus', 1.0078)
        self.mass_C_terminus = mass_config.get('mass_C_terminus', 17.0027)
        self.mass_CO = mass_config.get('mass_CO', 27.9949)
        self.mass_Phosphorylation = mass_config.get('mass_Phosphorylation', 79.96633)
        
        self.mass_AA = mass_config.get('mass_AA', {})
    
    def mAbdata_process(self, path: str, antibody: str) -> pd.DataFrame:
        """Process monoclonal antibody data."""
        # This would implement the mAbdata_process function from the notebook
        # For now, return empty DataFrame as placeholder
        return pd.DataFrame()
    
    def process_casanovo_v1(self, casanovo_path: str, mgf_file: str) -> pd.DataFrame:
        """Process CasanovoV1 results."""
        titles = []
        
        # Extract titles from MGF file
        try:
            with open(mgf_file, 'r') as fr:
                lines = fr.readlines()
                for line in lines:
                    if 'TITLE=' in line:
                        titles.append(line.strip().split('TITLE=')[-1])
        except FileNotFoundError:
            print(f"Warning: MGF file {mgf_file} not found.")
            return pd.DataFrame()
        
        # Find header line
        try:
            with open(casanovo_path) as f:
                for i, line in enumerate(f):
                    if line.startswith("PSH"):
                        casanovo_header = i
                        break
        except FileNotFoundError:
            print(f"Warning: Casanovo file {casanovo_path} not found.")
            return pd.DataFrame()
        
        # Load and process data
        df_deno = pd.read_csv(casanovo_path, sep="\t", header=casanovo_header)
        df_denovo = df_deno[df_deno['sequence'].notna()].copy()
        
        if df_denovo.empty:
            return pd.DataFrame()
        
        df_denovo['PSM_ID'] = df_denovo['PSM_ID'].astype(int)
        df_denovo['PSM_ID'] = df_denovo['PSM_ID'].apply(
            lambda x: titles[x] if x < len(titles) else None)
        
        # Process sequences and scores
        sequence = df_denovo['sequence'].tolist()
        PSM_ID = df_denovo['PSM_ID'].tolist()
        Score = [int((i * 100 + 100)/2) for i in df_denovo['search_engine_score[1]'].tolist()]
        
        # Process amino acid scores
        aascore = df_denovo['opt_ms_run[1]_aa_scores'].tolist()
        aascore = [str(i).replace(',', ' ') for i in aascore]
        aascore = [i.split() for i in aascore]
        aascore = [[str(int(float(j) * 100)) for j in i] for i in aascore]
        aascore = [" ".join(i) for i in aascore]
        
        return pd.DataFrame({
            'PSM_ID': PSM_ID,
            'sequence': sequence,
            'Score': Score,
            'aaScore': aascore
        })
    
    def process_casanovo_v2(self, casanovo_path: str, mgf_file: str) -> pd.DataFrame:
        """Process CasanovoV2 results."""
        # Similar implementation to CasanovoV1 but with different PSM_ID extraction
        titles = []
        
        # Extract titles from MGF file
        try:
            with open(mgf_file, 'r') as fr:
                lines = fr.readlines()
                for line in lines:
                    if 'TITLE=' in line:
                        titles.append(line.strip().split('TITLE=')[-1])
        except FileNotFoundError:
            print(f"Warning: MGF file {mgf_file} not found.")
            return pd.DataFrame()
        
        # Find header line
        try:
            with open(casanovo_path) as f:
                for i, line in enumerate(f):
                    if line.startswith("PSH"):
                        casanovo_header = i
                        break
        except FileNotFoundError:
            print(f"Warning: Casanovo file {casanovo_path} not found.")
            return pd.DataFrame()
        
        # Load and process data
        df_deno = pd.read_csv(casanovo_path, sep="\t", header=casanovo_header)
        df_denovo = df_deno[df_deno['sequence'].notna()].copy()
        
        if df_denovo.empty:
            return pd.DataFrame()
        
        df_denovo['PSM_ID'] = [int(item.replace('ms_run[1]:index=', '')) 
                              for item in df_denovo['spectra_ref']]
        df_denovo['PSM_ID'] = df_denovo['PSM_ID'].apply(
            lambda x: titles[x] if x < len(titles) else None)
        
        # Process sequences and scores (same as V1)
        sequence = df_denovo['sequence'].tolist()
        PSM_ID = df_denovo['PSM_ID'].tolist()
        Score = [int((i * 100 + 100)/2) for i in df_denovo['search_engine_score[1]'].tolist()]
        
        # Process amino acid scores
        aascore = df_denovo['opt_ms_run[1]_aa_scores'].tolist()
        aascore = [str(i).replace(',', ' ') for i in aascore]
        aascore = [i.split() for i in aascore]
        aascore = [[str(int(float(j) * 100)) for j in i] for i in aascore]
        aascore = [" ".join(i) for i in aascore]
        
        return pd.DataFrame({
            'PSM_ID': PSM_ID,
            'sequence': sequence,
            'Score': Score,
            'aaScore': aascore
        })
    
    def process_tool(self, tool: str, input_path: str, output_path: str, mgf_path: Optional[str] = None,
                    names_mapping: Optional[Dict[str, str]] = None) -> None:
        """
        Process a specific tool for all available samples.
        
        Args:
            tool: Tool name to process
            input_path: Base input path containing tool results
            output_path: Base output path for processed results
            mgf_path: Base path to MGF files
            names_mapping: Mapping from directory names to sample names
        """
        if tool not in self.config.tool_names:
            raise ValueError(f"Unknown tool: {tool}. Supported tools: {self.config.tool_names}")
        
        if names_mapping is None:
            names_mapping = self.config.default_names
        
        # Get tool-specific directory and file patterns from config
        tool_input_dirs = self.config.tool_input_dirs
        tool_file_patterns = self.config.tool_file_patterns
        
        tool_input_dir = tool_input_dirs.get(tool, tool.lower().replace('-', ''))
        tool_input_path = os.path.join(input_path, tool_input_dir)
        tool_output_path = os.path.join(output_path, tool)
        
        if not os.path.exists(tool_input_path):
            print(f"Warning: Input path {tool_input_path} does not exist")
            return
        
        os.makedirs(tool_output_path, exist_ok=True)
        mabs = os.listdir(tool_input_path)
        
        for mab in mabs:
            print(f"Processing {tool} for {mab}")
            mab_input_path = os.path.join(tool_input_path, mab)
            
            try:
                # Get file pattern for this tool
                file_pattern_info = tool_file_patterns.get(tool, {})
                input_file_pattern = file_pattern_info.get('input_file', '')
                input_subdir = file_pattern_info.get('input_subdir', '')
                requires_mgf = file_pattern_info.get('requires_mgf', False)
                
                # Build input file path
                if input_subdir:
                    input_file_base = os.path.join(mab_input_path, input_subdir)
                else:
                    input_file_base = mab_input_path
                
                # Format file pattern with variables
                sample_name = names_mapping.get(mab, mab)
                input_file_name = input_file_pattern.format(mab=mab, sample_name=sample_name)
                input_file = os.path.join(input_file_base, input_file_name)
                
                # Build MGF file path if needed
                mgf_file = None
                if requires_mgf and mgf_path:
                    mgf_pattern = self.config.mgf_file_pattern.format(sample_name=sample_name)
                    mgf_file = os.path.join(mgf_path, mab, self.config.mgf_subdir, mgf_pattern)
                
                # Call appropriate processing method
                result_df = pd.DataFrame()
                if tool == 'CasanovoV1':
                    result_df = self.process_casanovo_v1(input_file, mgf_file)
                elif tool == 'CasanovoV2':
                    result_df = self.process_casanovo_v2(input_file, mgf_file)
                # Add other tools as needed...
                else:
                    print(f"Tool {tool} processing not yet implemented")
                    continue
                
                # Save results
                if not result_df.empty:
                    output_file = os.path.join(tool_output_path, f"{mab}_processed.csv")
                    result_df.to_csv(output_file, index=False)
                    print(f"Saved results to {output_file}")
                
            except Exception as e:
                print(f"Error processing {tool} for {mab}: {e}")
                continue


def main():
    """Main function to run Denovo Process."""
    parser = argparse.ArgumentParser(description='Denovo Process Script')
    
    # Load config first to get default values
    config = get_denovo_config_loader()
    default_paths = config.default_paths
    
    parser.add_argument('--tool', type=str, required=True, choices=list(config.tool_names),
                       help='Tool to process')
    parser.add_argument('--input_path', type=str, required=True,
                       help='Base input path containing de novo results')
    parser.add_argument('--output_path', type=str,
                       default=default_paths.get('output_base', 'processed_results'),
                       help='Base output path for processed results')
    parser.add_argument('--mgf_path', type=str,
                       help='Base path to MGF files (required for some tools)')
    parser.add_argument('--names_mapping_file', type=str,
                       help='Path to JSON file containing names mapping')
    parser.add_argument('--config_path', type=str,
                       help='Path to configuration YAML file')

    args = parser.parse_args()

    # Load names mapping if provided
    names_mapping = None
    if args.names_mapping_file and os.path.exists(args.names_mapping_file):
        with open(args.names_mapping_file, 'r') as f:
            names_mapping = json.load(f)

    # Initialize Denovo processor with config
    processor = DenovoProcessor(config_path=args.config_path)

    # Create output directories
    os.makedirs(args.output_path, exist_ok=True)

    # Process the specified tool
    try:
        processor.process_tool(args.tool, args.input_path, args.output_path, 
                             args.mgf_path, names_mapping)
        print(f"Successfully processed {args.tool}")
    except Exception as e:
        print(f"Error processing {args.tool}: {e}")
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
