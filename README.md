# 🧬 **AbNovoBench: A Comprehensive, Standardized, and Reliable Benchmarking System for Evaluating Monoclonal Antibody De Novo Sequencing Analysis**

[![De Novo Sequencing](https://img.shields.io/badge/De%20Novo-Sequencing-blue)](https://github.com/dumbgoos/AbNovoBench)
[![Mass Spectrometry](https://img.shields.io/badge/Mass%20Spectrometry-Proteomics-green)](https://github.com/dumbgoos/AbNovoBench)
[![Python](https://img.shields.io/badge/Python-3.8%2B-brightgreen)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Website](https://img.shields.io/badge/Website-abnovobench.com-orange)](https://abnovobench.com)
[![Models](https://img.shields.io/badge/🤗%20Models-HuggingFace-yellow)](https://huggingface.co/LLMasterLL/AbNovobench)

---

## 🌟 **Introduction**

**AbNovoBench** is a comprehensive, standardized, and reliable benchmarking system for evaluating monoclonal antibody de novo sequencing analysis. This repository provides the first and largest high-quality dataset specifically designed for antibodies, along with standardized evaluation pipelines to compare multiple de novo sequencing tools.

### Key Features
- **Comprehensive Dataset**: The largest high-quality antibody MS dataset to date (1,638,248 PSMs from 131 mAbs)
- **Standardized Pipeline**: Automated processing and evaluation workflows
- **Multi-tool Support**: Evaluation of 13+ state-of-the-art de novo sequencing tools
- **Robust Metrics**: Multiple evaluation criteria including accuracy, robustness, and error analysis
- **Assembly Support**: Integration with ALPS and Stitch assembly algorithms
- **Online Platform**: Interactive benchmarking platform at [abnovobench.com](https://abnovobench.com)

### 🔗 **Important Links**

- **🌐 Official Website & Leaderboard**: [https://abnovobench.com](https://abnovobench.com)
  - Interactive platform with comprehensive benchmarking results
  - Real-time leaderboard of model performances  
  - Dataset exploration and visualization tools
  
- **🤗 Pre-trained Models**: [https://huggingface.co/LLMasterLL/AbNovobench](https://huggingface.co/LLMasterLL/AbNovobench)
  - Collection of 13 state-of-the-art de novo sequencing models
  - Ready-to-use pre-trained checkpoints
  - Model cards with detailed performance metrics

---

## 🚀 **Quick Start**

### Prerequisites
- Python 3.8 or higher
- Java Runtime Environment (for ALPS.jar)
- Required Python packages (see Installation)

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/dumbgoos/AbNovoBench.git
   cd AbNovoBench
   ```

2. Install Python dependencies using uv:
   ```bash
   # Install uv (if not already installed)
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # Install dependencies
   uv add pandas numpy PyYAML pyteomics tqdm scikit-learn matplotlib spectrum_utils fast-diff-match-patch
   ```

3. Download external tools:
   - **ALPS.jar**: Contact tool authors or obtain from their repository
   - **Stitch**: Download from [https://github.com/snijderlab/stitch](https://github.com/snijderlab/stitch)

---

## 📋 **Pipeline Overview**

AbNovoBench consists of a multi-stage pipeline for comprehensive evaluation:

```
1. Data Processing      →  2. Assembly (Optional)    →  3. Metric Analysis
   denovo_process.py       alps_fusion.py               accuracy_metric.py
   (Raw results → CSV)     assembly_stitch.py           robustness_metric.py
                          (CSV → Assembled seqs)        error_stats_metric.py
                                                        (Error pattern analysis)
```

### Stage 1: Data Processing ⭐ **REQUIRED FIRST STEP**
- **Script**: `src/summary/denovo_process.py`
- **Purpose**: Convert raw tool outputs to standardized CSV format
- **Input**: Raw de novo sequencing results from tools
- **Output**: 
  - Individual tool summary CSV files (`{Tool}_summary.csv`)
  - Merged summary file (`summary_merged.csv`) combining all tools
- **⚠️ Important**: This step MUST be completed before any metric analysis
- **Note**: The `summary_merged.csv` serves as the unified input for all downstream metric scripts

### Stage 2: Assembly (Optional)
- **Scripts**: `src/assembly/alps_fusion.py`, `src/assembly/assembly_stitch.py`
- **Purpose**: Improve sequence accuracy using assembly algorithms
- **Input**: Summary CSV files from Stage 1
- **Output**: Assembled sequences and improved predictions
- **Dependencies**: 
  - ALPS.jar executable (for alps_fusion.py)
  - Stitch binary (for assembly_stitch.py)

### Stage 3: Metric Analysis ⭐ **REQUIRES STAGE 1 OUTPUT**
- **Scripts**: `src/metric/accuracy_metric.py`, `src/metric/robustness_metric.py`, `src/metric/error_stats_metric.py`
- **Purpose**: Comprehensive evaluation and benchmarking
- **Input**: Summary CSV files from Stage 1 (and assembly results if available)
- **Output**: Detailed performance metrics and analysis reports
- **⚠️ Important**: These scripts require summary CSV files generated by `denovo_process.py`

---

## ⚡ **Execution Order - IMPORTANT**

**You MUST follow this order for the pipeline to work correctly:**

```bash
# Step 1: REQUIRED - Process raw results into summary CSV files
python src/summary/denovo_process.py --tool [TOOL_NAME] --input_path [RAW_RESULTS] --output_path [SUMMARY_OUTPUT]

# Step 2: OPTIONAL - Run assembly (requires summary files from Step 1)
python src/assembly/alps_fusion.py --input_path [SUMMARY_OUTPUT] ...
# OR
python src/assembly/assembly_stitch.py --input_path [SUMMARY_OUTPUT] ...

# Step 3: REQUIRED - Run metrics (requires summary files from Step 1)
python src/metric/accuracy_metric.py -i [SUMMARY_CSV] -o [METRICS_OUTPUT]
python src/metric/robustness_metric.py --csv-in [SUMMARY_CSV] --mgf-in [MGF_FILE] --csv-out [ROBUSTNESS_OUTPUT]
python src/metric/error_stats_metric.py --robustness_file [ROBUSTNESS_OUTPUT] --summary_file [SUMMARY_MERGED] --output_dir [ERROR_ANALYSIS]
```

**⚠️ Critical Notes:**
- `src/summary/denovo_process.py` MUST be run first to generate summary CSV files
- `src/metric/` scripts require the summary CSV files as input
- Assembly scripts also require summary CSV files, not raw tool outputs

**📁 File Dependencies:**
```
denovo_process.py generates:
├── {Tool}_summary.csv (for each tool)
└── summary_merged.csv (combined file for metrics)
    ├── Used by: accuracy_metric.py
    ├── Used by: robustness_metric.py  
    └── Used by: error_stats_metric.py (as --summary_file parameter)
```

---

## 🛠 **Detailed Usage Examples**

### 1. Process Raw Tool Results (FIRST STEP - REQUIRED)
```bash
python src/summary/denovo_process.py \
    --tool CasanovoV1 \
    --input_path /path/to/casanovo/raw/results \
    --output_path /path/to/summary/output
```

### 2. Run ALPS Assembly (OPTIONAL)
```bash
python src/assembly/alps_fusion.py \
    --tool CasanovoV1 \
    --input_path /path/to/summary/output \
    --output_path /path/to/alps/results \
    --mgf_path /path/to/mgf/files \
    --confidence_threshold_file data/Tool_Confidence_Threshold.csv
```

### 3. Run Stitch Assembly (OPTIONAL)
```bash
python src/assembly/assembly_stitch.py \
    --tool CasanovoV1 \
    --input_path /path/to/summary/output \
    --output_path /path/to/stitch/results \
    --mgf_path /path/to/mgf/files \
    --batchfiles_path /path/to/stitch/batchfiles
```

### 4. Calculate Accuracy Metrics
```bash
python src/metric/accuracy_metric.py \
    -i /path/to/summary.csv \
    -o /path/to/accuracy_results.csv
```

### 5. Calculate Robustness Metrics
```bash
python src/metric/robustness_metric.py \
    --mgf-in /path/to/spectra.mgf \
    --csv-in /path/to/summary.csv \
    --csv-out /path/to/robustness_results.csv
```

### 6. Analyze Error Statistics
```bash
python src/metric/error_stats_metric.py \
    --robustness_file /path/to/robustness_results.csv \
    --summary_file /path/to/summary_merged.csv \
    --output_dir /path/to/error_analysis
```
**Note**: The `--summary_file` parameter expects the `summary_merged.csv` file generated by `denovo_process.py` in Stage 1.

---

## 🛡 **Supported De Novo Sequencing Tools**

AbNovoBench currently supports evaluation of the following tools:

| Tool | Version | Reference |
|------|---------|-----------|
| **Casanovo** | v1, v2, v3 | [Noble-Lab/casanovo](https://github.com/Noble-Lab/casanovo) |
| **pi-HelixNovo** | Latest | [PHOENIXcenter/pi-HelixNovo](https://github.com/PHOENIXcenter/pi-HelixNovo) |
| **ContraNovo** | Latest | [BEAM-Labs/ContraNovo](https://github.com/BEAM-Labs/ContraNovo) |
| **InstaNovo** | Latest | [instadeepai/InstaNovo](https://github.com/instadeepai/InstaNovo) |
| **AdaNovo** | Latest | [Westlake-OmicsAI/adanovo_v1](https://github.com/Westlake-OmicsAI/adanovo_v1) |
| **DeepNovo** | Latest | [nh2tran/DeepNovo](https://github.com/nh2tran/DeepNovo) |
| **PointNovo** | Latest | [irleader/PointNovo](https://github.com/irleader/PointNovo) |
| **PGPointNovo** | Latest | [shallFun4Learning/PGPointNovo](https://github.com/shallFun4Learning/PGPointNovo) |
| **PepNet** | Latest | [lkytal/pepnet](https://github.com/lkytal/pepnet) |
| **NovoB** | Latest | [ProteomeTeam/NovoB](https://github.com/ProteomeTeam/NovoB) |
| **pNovo** | v3 | [pfind.org/software/pNovo](http://pfind.org/software/pNovo/index.html) |
| **SMSNet** | Latest | [cmb-chula/SMSNet](https://github.com/cmb-chula/SMSNet) |

---

## 📁 **Project Structure**

```
AbNovoBench/
├── src/
│   ├── assembly/           # Assembly algorithms
│   │   ├── alps_fusion.py      # ALPS assembly implementation
│   │   └── assembly_stitch.py  # Stitch assembly implementation
│   ├── config/            # Configuration files (YAML)
│   │   ├── alps_fusion.yaml
│   │   ├── assembly_stitch.yaml
│   │   ├── denovo_process.yaml
│   │   ├── error_stats_metric.yaml
│   │   ├── ms_noise_fragment_analysis.yaml
│   │   └── split_train_valid.yaml
│   ├── data/              # Data processing utilities
│   │   ├── ms_noise_fragment_analysis.py
│   │   └── split_train_valid.py
│   ├── jupyter/           # Original Jupyter notebooks (archived)
│   │   ├── alps_fusion.ipynb
│   │   ├── assembly_stitch.ipynb
│   │   ├── denovo_process.ipynb
│   │   ├── error_stats_metric.ipynb
│   │   ├── ms_noise_fragment_analysis.ipynb
│   │   └── split_train_valid.ipynb
│   ├── metric/            # Evaluation metrics
│   │   ├── accuracy_metric.py
│   │   ├── robustness_metric.py
│   │   └── error_stats_metric.py
│   └── summary/           # Data processing (START HERE)
│       └── denovo_process.py
├── envs/                 # Conda environment files for tools
└── README.md
```

---

## ⚙️ **Configuration**

All scripts use YAML configuration files located in `src/config/`. These files contain:
- Tool-specific parameters
- File path configurations  
- Processing parameters
- Default values

Modify the configuration files to adapt the pipeline to your specific setup and data organization.

---

## 📊 **Evaluation Metrics**

### Accuracy Metrics
- **Amino Acid Precision/Recall**: Character-level accuracy
- **Peptide-level Accuracy**: Exact sequence match rates
- **PTM Accuracy**: Post-translational modification detection
- **AUC Scores**: Ranking quality assessment

### Robustness Metrics  
- **Fragment Ion Coverage**: Spectral explanation quality
- **Noise Factor Analysis**: Signal-to-noise ratio evaluation
- **Spectral Quality Scores**: MS/MS data quality assessment

### Error Analysis
- **Substitution Errors**: Amino acid replacement patterns
- **Insertion/Deletion Errors**: Sequence length variations
- **Permutation Errors**: Sequence rearrangement detection

---

## 🔧 **Dependencies and External Tools**

### Required External Tools
1. **ALPS.jar**: Assembly algorithm (contact authors for access)
2. **Stitch**: Template-based assembly tool
   - Download: [https://github.com/snijderlab/stitch](https://github.com/snijderlab/stitch)
   - Distributed executables available for Windows, Linux, and macOS

### Pre-trained Models
- **Download all models**: [https://huggingface.co/LLMasterLL/AbNovobench](https://huggingface.co/LLMasterLL/AbNovobench)
- **Supported tools**: AdaNovo, CasaNovo (V1/V2), ContraNovo, DeepNovo, InstaNovo, PepNet, PGPointNovo, pi-HelixNovo, pi-PrimeNovo, PointNovo, SMSNet
- **Model formats**: PyTorch checkpoints (.ckpt, .pth), TensorFlow models (.h5)

### Online Resources
- **Official Website**: [https://abnovobench.com](https://abnovobench.com)
  - Interactive leaderboard with real-time performance comparisons
  - Dataset exploration and visualization tools
  - Comprehensive benchmarking results and analysis
  - Submit your own models for evaluation

### Python Dependencies
- pandas, numpy: Data manipulation
- PyYAML: Configuration file parsing
- pyteomics: Mass spectrometry data handling
- scikit-learn: Machine learning metrics
- matplotlib: Visualization
- spectrum_utils: Spectral analysis utilities
- fast_diff_match_patch: Text difference analysis (for error_stats_metric)

---

## 🤝 **Contributing**

We welcome contributions to improve AbNovoBench! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📚 **Citation**

If you use AbNovoBench in your research, please cite:

```bibtex
@misc{jiang2025abnovobench,
  title        = {AbNovoBench: A Comprehensive, Standardized, and Reliable Benchmarking System for Evaluating Monoclonal Antibody De Novo Sequencing Analysis},
  author       = {Wenbin Jiang and Ling Luo and Lihong Huang and Jin Xiao and Zihan Lin and Yijie Qiu and Jiying Wang and Ouyang Hu and Sainan Zhang and Mengsha Tong and Ningshao Xia and Yueting Xiong and Quan Yuan and Rongshan Yu},
  year         = {2025},
  howpublished = {https://github.com/dumbgoos/AbNovoBench}
}
```

---

## 📜 **License**

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

---

## 📞 **Support**

For questions, issues, or support:
- Open an issue on [GitHub Issues](https://github.com/dumbgoos/AbNovoBench/issues)
- Contact the authors via the repository

---

**AbNovoBench** - Advancing antibody de novo sequencing through comprehensive benchmarking and standardized evaluation.