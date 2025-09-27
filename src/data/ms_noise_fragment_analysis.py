#!/usr/bin/env python3
"""
MS Noise Fragment Analysis Script

This script processes mass spectrometry (MS) data to calculate fragment ion information
and noise factors. It uses KDTree-based matching for efficient m/z value matching
and supports multiprocessing for parallel processing of spectra.

Usage:
    python MS_Noise_Fragment_Analysis.py --input_file /path/to/input.mgf --output_dir /path/to/output
"""

import os
import argparse
import re
import yaml
import json
import numpy as np
import pandas as pd
from typing import Dict, Optional, Any, List, Tuple
from pyteomics import mgf
from tqdm import tqdm
from multiprocessing import Pool
from scipy.spatial import cKDTree


class MSAnalysisConfigLoader:
    """Configuration loader for MS Noise Fragment Analysis parameters."""
    
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
            config_path = os.path.join(src_dir, 'config', 'ms_noise_fragment_analysis.yaml')
        
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
            'analysis_parameters': {
                'mass_tolerance': 0.02,
                'use_kdtree': True,
                'parallel_processing': True,
                'show_progress': True
            },
            'multiprocessing_config': {
                'num_processes': 4,
                'chunk_size': 100
            },
            'default_paths': {
                'input_base': 'data',
                'output_base': 'analysis_results'
            }
        }
    
    @property
    def analysis_parameters(self) -> Dict[str, Any]:
        """Get analysis parameters."""
        return self.config.get('analysis_parameters', {})
    
    @property
    def multiprocessing_config(self) -> Dict[str, Any]:
        """Get multiprocessing configuration."""
        return self.config.get('multiprocessing_config', {})
    
    @property
    def progress_config(self) -> Dict[str, Any]:
        """Get progress tracking configuration."""
        return self.config.get('progress_config', {})
    
    @property
    def fragment_ion_types(self) -> Dict[str, bool]:
        """Get fragment ion types configuration."""
        return self.config.get('fragment_ion_types', {})
    
    @property
    def mass_constants(self) -> Dict[str, Any]:
        """Get mass constants."""
        return self.config.get('mass_constants', {})
    
    @property
    def modification_masses(self) -> Dict[str, float]:
        """Get modification masses."""
        return self.config.get('modification_masses', {})
    
    @property
    def noise_analysis(self) -> Dict[str, Any]:
        """Get noise analysis parameters."""
        return self.config.get('noise_analysis', {})
    
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
    def output_config(self) -> Dict[str, Any]:
        """Get output configuration."""
        return self.config.get('output_config', {})
    
    @property
    def kdtree_config(self) -> Dict[str, Any]:
        """Get KDTree configuration."""
        return self.config.get('kdtree_config', {})
    
    @property
    def sequence_processing(self) -> Dict[str, Any]:
        """Get sequence processing configuration."""
        return self.config.get('sequence_processing', {})
    
    @property
    def quality_filters(self) -> Dict[str, Any]:
        """Get quality control filters."""
        return self.config.get('quality_filters', {})


# Global configuration loader
_ms_analysis_config_loader = None

def get_ms_analysis_config_loader(config_path: Optional[str] = None) -> MSAnalysisConfigLoader:
    """Get global configuration loader instance."""
    global _ms_analysis_config_loader
    if _ms_analysis_config_loader is None or config_path is not None:
        _ms_analysis_config_loader = MSAnalysisConfigLoader(config_path)
    return _ms_analysis_config_loader


class MSNoiseFragmentAnalyzer:
    """MS noise and fragment ion analyzer."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the analyzer.
        
        Args:
            config_path: Path to configuration file
        """
        # Load configuration
        self.config = get_ms_analysis_config_loader(config_path)
        
        # Set parameters from config
        analysis_params = self.config.analysis_parameters
        self.mass_tolerance = analysis_params.get('mass_tolerance', 0.02)
        self.use_kdtree = analysis_params.get('use_kdtree', True)
        self.parallel_processing = analysis_params.get('parallel_processing', True)
        self.show_progress = analysis_params.get('show_progress', True)
        
        # Multiprocessing config
        mp_config = self.config.multiprocessing_config
        self.num_processes = mp_config.get('num_processes', 4)
        self.chunk_size = mp_config.get('chunk_size', 100)
        
        # Load mass constants
        mass_constants = self.config.mass_constants
        self.proton_mass = mass_constants.get('proton_mass', 1.00727647)
        self.water_mass = mass_constants.get('water_mass', 18.0105647)
        self.ammonia_mass = mass_constants.get('ammonia_mass', 17.0265491)
        self.amino_acid_masses = mass_constants.get('amino_acid_masses', {})
        
        # Load modification masses
        self.modification_masses = self.config.modification_masses
        
        # Fragment ion types
        self.fragment_ion_types = self.config.fragment_ion_types
    
    def extract_modified_amino_acids(self, peptide_sequence: str) -> List[Tuple[int, str]]:
        """
        Extracts modified amino acids from a peptide sequence.
        
        Args:
            peptide_sequence: Peptide sequence with modifications
            
        Returns:
            List of tuples containing (position, modified_amino_acid)
        """
        matches = re.finditer(r'[A-Za-z]\([^)]+\)', peptide_sequence)
        results = [(match.start(), match.group()) for match in matches]
        return results
    
    def calculate_theoretical_fragments(self, peptide_sequence: str, charge: int = 1) -> Dict[str, List[float]]:
        """
        Calculate theoretical fragment ion m/z values for a peptide.
        
        Args:
            peptide_sequence: Peptide sequence (with modifications)
            charge: Precursor charge state
            
        Returns:
            Dictionary containing fragment ion types and their m/z values
        """
        # Parse modifications
        modifications = self.extract_modified_amino_acids(peptide_sequence)
        
        # Clean sequence (remove modification annotations)
        clean_sequence = re.sub(r'\([^)]+\)', '', peptide_sequence)
        
        # Calculate amino acid masses including modifications
        aa_masses = []
        mod_positions = {pos: mod for pos, mod in modifications}
        
        current_pos = 0
        for i, aa in enumerate(clean_sequence):
            mass = self.amino_acid_masses.get(aa, 0.0)
            
            # Check for modifications at this position
            for mod_pos, mod_aa in mod_positions.items():
                if mod_pos <= current_pos:
                    mod_mass = self.modification_masses.get(mod_aa, 0.0)
                    mass += mod_mass
                    break
            
            aa_masses.append(mass)
            current_pos += len(aa)
        
        # Calculate fragment ions
        fragments = {}
        
        if self.fragment_ion_types.get('b_ions', True):
            fragments['b_ions'] = self._calculate_b_ions(aa_masses, charge)
        
        if self.fragment_ion_types.get('y_ions', True):
            fragments['y_ions'] = self._calculate_y_ions(aa_masses, charge)
        
        if self.fragment_ion_types.get('a_ions', False):
            fragments['a_ions'] = self._calculate_a_ions(aa_masses, charge)
        
        return fragments
    
    def _calculate_b_ions(self, aa_masses: List[float], charge: int) -> List[float]:
        """Calculate b-ion m/z values."""
        b_ions = []
        cumulative_mass = 0.0
        
        for i in range(len(aa_masses) - 1):  # b-ions go up to n-1
            cumulative_mass += aa_masses[i]
            # b-ion = cumulative mass + proton mass
            mz = (cumulative_mass + self.proton_mass * charge) / charge
            b_ions.append(mz)
        
        return b_ions
    
    def _calculate_y_ions(self, aa_masses: List[float], charge: int) -> List[float]:
        """Calculate y-ion m/z values."""
        y_ions = []
        cumulative_mass = self.water_mass  # y-ions include water
        
        for i in range(len(aa_masses) - 1, 0, -1):  # y-ions go from n-1 to 1
            cumulative_mass += aa_masses[i]
            # y-ion = cumulative mass + proton mass
            mz = (cumulative_mass + self.proton_mass * charge) / charge
            y_ions.append(mz)
        
        return y_ions
    
    def _calculate_a_ions(self, aa_masses: List[float], charge: int) -> List[float]:
        """Calculate a-ion m/z values."""
        a_ions = []
        cumulative_mass = 0.0
        
        for i in range(len(aa_masses) - 1):  # a-ions go up to n-1
            cumulative_mass += aa_masses[i]
            # a-ion = b-ion - CO
            mz = (cumulative_mass + self.proton_mass * charge - 27.9949146) / charge
            a_ions.append(mz)
        
        return a_ions
    
    def match_fragments_kdtree(self, observed_mz: np.ndarray, theoretical_mz: List[float]) -> Dict[str, Any]:
        """
        Match observed m/z values to theoretical fragments using KDTree.
        
        Args:
            observed_mz: Array of observed m/z values
            theoretical_mz: List of theoretical m/z values
            
        Returns:
            Dictionary containing matching statistics
        """
        if not theoretical_mz:
            return {'matched_count': 0, 'total_theoretical': 0, 'match_percentage': 0.0}
        
        # Build KDTree for theoretical m/z values
        theoretical_array = np.array(theoretical_mz).reshape(-1, 1)
        kdtree_config = self.config.kdtree_config
        leaf_size = kdtree_config.get('leaf_size', 30)
        
        tree = cKDTree(theoretical_array, leafsize=leaf_size)
        
        # Query for matches within tolerance
        matches = 0
        for mz in observed_mz:
            distances, indices = tree.query([[mz]], k=1, distance_upper_bound=self.mass_tolerance)
            if distances[0] < np.inf:  # Found a match within tolerance
                matches += 1
        
        match_percentage = (matches / len(theoretical_mz)) * 100 if theoretical_mz else 0.0
        
        return {
            'matched_count': matches,
            'total_theoretical': len(theoretical_mz),
            'match_percentage': match_percentage
        }
    
    def analyze_noise(self, spectrum: Dict[str, Any]) -> Dict[str, float]:
        """
        Analyze noise characteristics of a spectrum.
        
        Args:
            spectrum: Spectrum dictionary from pyteomics
            
        Returns:
            Dictionary containing noise statistics
        """
        intensities = spectrum['intensity array']
        
        noise_params = self.config.noise_analysis
        intensity_threshold = noise_params.get('intensity_threshold', 0.01)
        noise_percentile = noise_params.get('noise_percentile', 10)
        signal_to_noise_ratio = noise_params.get('signal_to_noise_ratio', 3.0)
        
        # Normalize intensities
        max_intensity = np.max(intensities)
        normalized_intensities = intensities / max_intensity
        
        # Calculate noise floor
        noise_floor = np.percentile(normalized_intensities, noise_percentile)
        
        # Calculate signal-to-noise ratio
        signal_peaks = normalized_intensities[normalized_intensities > intensity_threshold]
        avg_signal = np.mean(signal_peaks) if len(signal_peaks) > 0 else 0.0
        snr = avg_signal / noise_floor if noise_floor > 0 else 0.0
        
        return {
            'noise_floor': noise_floor,
            'signal_to_noise_ratio': snr,
            'num_signal_peaks': len(signal_peaks),
            'avg_signal_intensity': avg_signal
        }
    
    def process_spectrum(self, spectrum: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single spectrum for fragment analysis.
        
        Args:
            spectrum: Spectrum dictionary from pyteomics
            
        Returns:
            Dictionary containing analysis results
        """
        try:
            # Extract basic information
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']
            params = spectrum.get('params', {})
            
            # Apply quality filters
            quality_filters = self.config.quality_filters
            min_peaks = quality_filters.get('min_peaks', 10)
            min_intensity = quality_filters.get('min_intensity', 0.001)
            
            if len(mz_array) < min_peaks:
                return {'error': 'Too few peaks'}
            
            if np.max(intensity_array) < min_intensity:
                return {'error': 'Intensity too low'}
            
            # Get peptide sequence
            peptide_seq = params.get('seq', '')
            if not peptide_seq:
                return {'error': 'No peptide sequence'}
            
            # Get charge state
            charge = params.get('charge', [1])[0] if isinstance(params.get('charge', [1]), list) else params.get('charge', 1)
            
            # Calculate theoretical fragments
            theoretical_fragments = self.calculate_theoretical_fragments(peptide_seq, charge)
            
            # Match fragments
            results = {'peptide_sequence': peptide_seq, 'charge': charge}
            
            for fragment_type, theoretical_mz in theoretical_fragments.items():
                match_results = self.match_fragments_kdtree(mz_array, theoretical_mz)
                results[fragment_type] = match_results
            
            # Analyze noise
            noise_results = self.analyze_noise(spectrum)
            results['noise_analysis'] = noise_results
            
            return results
            
        except Exception as e:
            return {'error': str(e)}
    
    def analyze_mgf_file(self, input_file: str, output_dir: str) -> str:
        """
        Analyze an MGF file for fragment ions and noise.
        
        Args:
            input_file: Path to input MGF file
            output_dir: Output directory for results
            
        Returns:
            Path to output results file
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Determine output file name
        base_name = os.path.basename(input_file).replace('.mgf', '')
        output_config = self.config.output_config
        result_suffix = output_config.get('result_suffix', '_fragment_analysis')
        output_format = output_config.get('output_format', 'csv')
        output_file = os.path.join(output_dir, f"{base_name}{result_suffix}.{output_format}")
        
        # Read spectra
        file_params = self.config.file_parameters
        with open(input_file, 'r') as f:
            spectra = mgf.read(f, 
                             convert_arrays=file_params.get('convert_arrays', True),
                             read_charges=file_params.get('read_charges', True),
                             dtype=file_params.get('dtype', 'float32'),
                             use_index=file_params.get('use_index', False))
            
            spectra_list = list(spectra)
        
        print(f"Processing {len(spectra_list)} spectra from {input_file}")
        
        # Process spectra
        if self.parallel_processing and len(spectra_list) > self.chunk_size:
            with Pool(processes=self.num_processes) as pool:
                if self.show_progress:
                    results = list(tqdm(
                        pool.imap(self.process_spectrum, spectra_list, chunksize=self.chunk_size),
                        total=len(spectra_list),
                        desc="Processing spectra"
                    ))
                else:
                    results = pool.map(self.process_spectrum, spectra_list)
        else:
            if self.show_progress:
                results = [self.process_spectrum(spectrum) for spectrum in tqdm(spectra_list, desc="Processing spectra")]
            else:
                results = [self.process_spectrum(spectrum) for spectrum in spectra_list]
        
        # Save results
        results_df = pd.DataFrame(results)
        
        if output_format.lower() == 'csv':
            results_df.to_csv(output_file, index=False)
        elif output_format.lower() == 'tsv':
            results_df.to_csv(output_file, sep='\t', index=False)
        elif output_format.lower() == 'json':
            results_df.to_json(output_file, orient='records', indent=2)
        
        print(f"Results saved to {output_file}")
        return output_file
    
    def batch_analyze(self, input_dir: str, output_dir: str, file_pattern: str = "*.mgf") -> None:
        """
        Analyze multiple MGF files in a directory.
        
        Args:
            input_dir: Input directory containing MGF files
            output_dir: Output directory for results
            file_pattern: File pattern to match (default: "*.mgf")
        """
        import glob
        
        if not os.path.exists(input_dir):
            raise FileNotFoundError(f"Input directory {input_dir} not found")
        
        # Find all MGF files
        search_pattern = os.path.join(input_dir, file_pattern)
        mgf_files = glob.glob(search_pattern)
        
        if not mgf_files:
            print(f"No MGF files found in {input_dir} matching pattern {file_pattern}")
            return
        
        print(f"Found {len(mgf_files)} MGF files to analyze")
        
        for mgf_file in mgf_files:
            print(f"\nAnalyzing {os.path.basename(mgf_file)}...")
            try:
                self.analyze_mgf_file(mgf_file, output_dir)
            except Exception as e:
                print(f"Error analyzing {mgf_file}: {e}")


def main():
    """Main function to run MS Noise Fragment Analysis."""
    parser = argparse.ArgumentParser(description='MS Noise Fragment Analysis Script')
    
    # Load config first to get default values
    config = get_ms_analysis_config_loader()
    analysis_params = config.analysis_parameters
    default_paths = config.default_paths
    
    parser.add_argument('--input_file', type=str,
                       help='Input MGF file to analyze')
    parser.add_argument('--input_dir', type=str,
                       help='Input directory containing MGF files for batch analysis')
    parser.add_argument('--output_dir', type=str,
                       default=default_paths.get('output_base', 'analysis_results'),
                       help='Output directory for analysis results')
    parser.add_argument('--mass_tolerance', type=float,
                       default=analysis_params.get('mass_tolerance', 0.02),
                       help='Mass tolerance for fragment matching (Da)')
    parser.add_argument('--num_processes', type=int,
                       default=config.multiprocessing_config.get('num_processes', 4),
                       help='Number of parallel processes')
    parser.add_argument('--disable_parallel', action='store_true',
                       help='Disable parallel processing')
    parser.add_argument('--config_path', type=str,
                       help='Path to configuration YAML file')
    parser.add_argument('--file_pattern', type=str, default="*.mgf",
                       help='File pattern for batch processing (default: *.mgf)')

    args = parser.parse_args()

    if not args.input_file and not args.input_dir:
        parser.error("Either --input_file or --input_dir must be specified")

    # Initialize analyzer with config
    analyzer = MSNoiseFragmentAnalyzer(config_path=args.config_path)
    
    # Override parameters if provided
    analyzer.mass_tolerance = args.mass_tolerance
    analyzer.num_processes = args.num_processes
    if args.disable_parallel:
        analyzer.parallel_processing = False

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        if args.input_file:
            # Process single file
            output_file = analyzer.analyze_mgf_file(args.input_file, args.output_dir)
            print(f"Successfully analyzed {args.input_file}")
        elif args.input_dir:
            # Process directory
            analyzer.batch_analyze(args.input_dir, args.output_dir, args.file_pattern)
            print(f"Successfully analyzed all files in {args.input_dir}")
            
    except Exception as e:
        print(f"Error during analysis: {e}")
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
