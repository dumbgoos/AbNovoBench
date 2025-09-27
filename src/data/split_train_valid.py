#!/usr/bin/env python3
"""
Split Train Valid Script

This script splits MGF files into training and validation subsets. It ensures that
peptides in the validation set are unique and do not overlap with the training set.

Usage:
    python Split_Train_Valid.py --input_file /path/to/input.mgf --validation_ratio 0.1
"""

import os
import argparse
import random
import yaml
import json
from typing import Dict, Optional, Any, List
from pyteomics import mgf
from tqdm import tqdm


class SplitConfigLoader:
    """Configuration loader for Split Train Valid parameters."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration loader.
        
        Args:
            config_path: Path to configuration file. If None, uses default path.
        """
        if config_path is None:
            # Default config path is in the src/config directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            # Go up one level from data to src, then into config
            src_dir = os.path.dirname(script_dir)
            config_path = os.path.join(src_dir, 'config', 'split_train_valid.yaml')
        
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
        # Go up two levels: data -> src -> project_root
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
            'split_parameters': {
                'validation_ratio': 0.1,
                'random_seed': 42,
                'ensure_unique_sequences': True
            },
            'file_parameters': {
                'convert_arrays': True,
                'read_charges': False,
                'dtype': 'float32',
                'use_index': False
            },
            'default_paths': {
                'input_base': 'data',
                'output_base': 'data/split'
            },
            'file_patterns': {
                'input_extension': '.mgf',
                'train_suffix': '_train.mgf',
                'valid_suffix': '_valid.mgf'
            }
        }
    
    @property
    def split_parameters(self) -> Dict[str, Any]:
        """Get split parameters."""
        return self.config.get('split_parameters', {})
    
    @property
    def file_parameters(self) -> Dict[str, Any]:
        """Get file processing parameters."""
        return self.config.get('file_parameters', {})
    
    @property
    def default_paths(self) -> Dict[str, str]:
        """Get default file paths (resolved to absolute paths)."""
        paths = self.config.get('default_paths', {})
        resolved_paths = {}
        for key, path in paths.items():
            resolved_paths[key] = self._resolve_path(path)
        return resolved_paths
    
    @property
    def file_patterns(self) -> Dict[str, str]:
        """Get file patterns."""
        return self.config.get('file_patterns', {})
    
    @property
    def progress_config(self) -> Dict[str, Any]:
        """Get progress tracking configuration."""
        return self.config.get('progress_config', {})
    
    @property
    def sequence_filters(self) -> Dict[str, Any]:
        """Get sequence filtering parameters."""
        return self.config.get('sequence_filters', {})
    
    @property
    def output_format(self) -> Dict[str, Any]:
        """Get output format configuration."""
        return self.config.get('output_format', {})


# Global configuration loader
_split_config_loader = None

def get_split_config_loader(config_path: Optional[str] = None) -> SplitConfigLoader:
    """Get global configuration loader instance."""
    global _split_config_loader
    if _split_config_loader is None or config_path is not None:
        _split_config_loader = SplitConfigLoader(config_path)
    return _split_config_loader


class TrainValidSplitter:
    """Train/validation data splitter for MGF files."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the splitter.
        
        Args:
            config_path: Path to configuration file
        """
        # Load configuration
        self.config = get_split_config_loader(config_path)
        
        # Set parameters from config
        split_params = self.config.split_parameters
        self.validation_ratio = split_params.get('validation_ratio', 0.1)
        self.random_seed = split_params.get('random_seed', 42)
        self.ensure_unique_sequences = split_params.get('ensure_unique_sequences', True)
        
        file_params = self.config.file_parameters
        self.convert_arrays = file_params.get('convert_arrays', True)
        self.read_charges = file_params.get('read_charges', False)
        self.dtype = file_params.get('dtype', 'float32')
        self.use_index = file_params.get('use_index', False)
        
        # Set random seed for reproducibility
        random.seed(self.random_seed)
    
    def split_for_training(self, input_file: str, output_dir: Optional[str] = None) -> tuple:
        """
        Splits an MGF file into training and validation subsets.
        
        Args:
            input_file: The input MGF file path
            output_dir: Output directory (if None, uses input file directory)
            
        Returns:
            tuple: (train_file_path, valid_file_path)
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found")
        
        # Determine output paths
        if output_dir is None:
            output_dir = os.path.dirname(input_file)
        else:
            os.makedirs(output_dir, exist_ok=True)
        
        file_patterns = self.config.file_patterns
        base_name = os.path.basename(input_file).replace(file_patterns.get('input_extension', '.mgf'), '')
        train_file = os.path.join(output_dir, base_name + file_patterns.get('train_suffix', '_train.mgf'))
        valid_file = os.path.join(output_dir, base_name + file_patterns.get('valid_suffix', '_valid.mgf'))
        
        # Open files for reading and writing
        with open(input_file, "r") as file1, \
             open(valid_file, "w") as valid_f, \
             open(train_file, "w") as training_f:
            
            # Read spectra from MGF file
            progress_config = self.config.progress_config
            show_progress = progress_config.get('show_progress', True)
            progress_desc = progress_config.get('progress_desc', 'Reading spectra')
            progress_unit = progress_config.get('progress_unit', 'spectrum')
            
            sps = mgf.read(file1, 
                          convert_arrays=self.convert_arrays,
                          read_charges=self.read_charges,
                          dtype=self.dtype,
                          use_index=self.use_index)
            
            if show_progress:
                list_of_spectra = [sp for sp in tqdm(sps, desc=progress_desc, unit=progress_unit)]
            else:
                list_of_spectra = list(sps)
            
            # Shuffle spectra to randomize distribution
            random.shuffle(list_of_spectra)
            
            # Extract sequences and apply filters
            filtered_spectra = self._apply_sequence_filters(list_of_spectra)
            
            if self.ensure_unique_sequences:
                train_spectra, valid_spectra = self._split_with_unique_sequences(filtered_spectra)
            else:
                train_spectra, valid_spectra = self._split_simple(filtered_spectra)
            
            # Write training spectra
            if show_progress:
                train_iter = tqdm(train_spectra, desc="Writing training spectra", unit="spectrum")
            else:
                train_iter = train_spectra
            
            for sp in train_iter:
                mgf.write([sp], training_f)
            
            # Write validation spectra
            if show_progress:
                valid_iter = tqdm(valid_spectra, desc="Writing validation spectra", unit="spectrum")
            else:
                valid_iter = valid_spectra
            
            for sp in valid_iter:
                mgf.write([sp], valid_f)
        
        print(f"Split completed:")
        print(f"  Training set: {len(train_spectra)} spectra -> {train_file}")
        print(f"  Validation set: {len(valid_spectra)} spectra -> {valid_file}")
        
        return train_file, valid_file
    
    def _apply_sequence_filters(self, spectra: List[Dict]) -> List[Dict]:
        """Apply sequence filtering based on configuration."""
        sequence_filters = self.config.sequence_filters
        
        if not sequence_filters:
            return spectra
        
        min_length = sequence_filters.get('min_length', 0)
        max_length = sequence_filters.get('max_length', float('inf'))
        exclude_modifications = sequence_filters.get('exclude_modifications', False)
        
        filtered_spectra = []
        
        for sp in spectra:
            if 'seq' not in sp['params']:
                continue
            
            sequence = sp['params']['seq']
            
            # Check sequence length
            clean_seq = sequence.replace('(', '').replace(')', '')  # Remove modification markers for length check
            if len(clean_seq) < min_length or len(clean_seq) > max_length:
                continue
            
            # Check modifications
            if exclude_modifications and ('(' in sequence or ')' in sequence):
                continue
            
            filtered_spectra.append(sp)
        
        return filtered_spectra
    
    def _split_with_unique_sequences(self, spectra: List[Dict]) -> tuple:
        """Split ensuring unique sequences in validation set."""
        # Extract unique sequences from spectra
        seq_all = [sp['params']['seq'] for sp in spectra if 'seq' in sp['params']]
        seq_unique = list(set(seq_all))
        
        # Calculate number of unique sequences for validation
        num_valid_sequences = int(len(seq_unique) * self.validation_ratio)
        
        # Select unique sequences for validation
        random.shuffle(seq_unique)
        valid_sequences = set(seq_unique[:num_valid_sequences])
        
        # Split spectra based on sequence membership
        train_spectra = []
        valid_spectra = []
        
        for sp in spectra:
            if 'seq' in sp['params'] and sp['params']['seq'] in valid_sequences:
                valid_spectra.append(sp)
            else:
                train_spectra.append(sp)
        
        return train_spectra, valid_spectra
    
    def _split_simple(self, spectra: List[Dict]) -> tuple:
        """Simple split without ensuring unique sequences."""
        num_valid = int(len(spectra) * self.validation_ratio)
        
        valid_spectra = spectra[:num_valid]
        train_spectra = spectra[num_valid:]
        
        return train_spectra, valid_spectra
    
    def batch_split(self, input_dir: str, output_dir: str, file_pattern: str = "*.mgf") -> None:
        """
        Split multiple MGF files in a directory.
        
        Args:
            input_dir: Input directory containing MGF files
            output_dir: Output directory for split files
            file_pattern: File pattern to match (default: "*.mgf")
        """
        import glob
        
        if not os.path.exists(input_dir):
            raise FileNotFoundError(f"Input directory {input_dir} not found")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Find all MGF files
        search_pattern = os.path.join(input_dir, file_pattern)
        mgf_files = glob.glob(search_pattern)
        
        if not mgf_files:
            print(f"No MGF files found in {input_dir} matching pattern {file_pattern}")
            return
        
        print(f"Found {len(mgf_files)} MGF files to process")
        
        for mgf_file in mgf_files:
            print(f"\nProcessing {os.path.basename(mgf_file)}...")
            try:
                self.split_for_training(mgf_file, output_dir)
            except Exception as e:
                print(f"Error processing {mgf_file}: {e}")


def main():
    """Main function to run Train/Valid split."""
    parser = argparse.ArgumentParser(description='Split Train Valid Script')
    
    # Load config first to get default values
    config = get_split_config_loader()
    split_params = config.split_parameters
    default_paths = config.default_paths
    
    parser.add_argument('--input_file', type=str,
                       help='Input MGF file to split')
    parser.add_argument('--input_dir', type=str,
                       help='Input directory containing MGF files for batch processing')
    parser.add_argument('--output_dir', type=str,
                       default=default_paths.get('output_base', 'data/split'),
                       help='Output directory for split files')
    parser.add_argument('--validation_ratio', type=float,
                       default=split_params.get('validation_ratio', 0.1),
                       help='Ratio of data for validation (default: 0.1)')
    parser.add_argument('--random_seed', type=int,
                       default=split_params.get('random_seed', 42),
                       help='Random seed for reproducible splits')
    parser.add_argument('--ensure_unique_sequences', action='store_true',
                       default=split_params.get('ensure_unique_sequences', True),
                       help='Ensure validation sequences are unique')
    parser.add_argument('--config_path', type=str,
                       help='Path to configuration YAML file')
    parser.add_argument('--file_pattern', type=str, default="*.mgf",
                       help='File pattern for batch processing (default: *.mgf)')

    args = parser.parse_args()

    if not args.input_file and not args.input_dir:
        parser.error("Either --input_file or --input_dir must be specified")

    # Initialize splitter with config
    splitter = TrainValidSplitter(config_path=args.config_path)
    
    # Override parameters if provided
    splitter.validation_ratio = args.validation_ratio
    splitter.random_seed = args.random_seed
    splitter.ensure_unique_sequences = args.ensure_unique_sequences
    
    # Set random seed
    random.seed(splitter.random_seed)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        if args.input_file:
            # Process single file
            train_file, valid_file = splitter.split_for_training(args.input_file, args.output_dir)
            print(f"Successfully split {args.input_file}")
        elif args.input_dir:
            # Process directory
            splitter.batch_split(args.input_dir, args.output_dir, args.file_pattern)
            print(f"Successfully processed all files in {args.input_dir}")
            
    except Exception as e:
        print(f"Error during splitting: {e}")
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
