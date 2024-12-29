# 🧬 **Antibody Mass Spectrometry De Novo Benchmark**

Welcome to the **Antibody Mass Spectrometry De Novo Benchmark** repository! This project aims to evaluate and benchmark de novo sequencing tools specifically for antibody mass spectrometry (MS) data analysis. 

![De Novo Sequencing](https://img.shields.io/badge/De%20Novo-Sequencing-blue)
![Mass Spectrometry](https://img.shields.io/badge/Mass%20Spectrometry-Proteomics-green)
---

## 📝 **TODO**

- [x] Upload the environment yml file of the supported denovo method
- [ ] Release Docker
- [ ] Release benchmark pipeline

---

## 🌟 **Introduction**

Mass spectrometry-based de novo sequencing is a critical method for understanding antibody structures and functionalities. This repository provides:

- A curated dataset for benchmarking.
- Scripts and pipelines to compare multiple de novo sequencing tools.
- Visualized results for evaluating tool performance.

---

## 🛠 **Installation**

1. Clone the repository:
   ```bash
   git clone https://github.com/dumbgoos/AbNovoBench.git
   ```

2. Navigate to the repository:
   ```bash
   cd AbNovoBench
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---

## 🚀 **Usage**

1. Prepare your MS data in the required format.
2. Run the benchmarking pipeline:
   ```bash
   python benchmark.py --input /path/to/data --tools tool1,tool2 --output /path/to/results
   ```
3. Visualize results:
   ```bash
   python visualize_results.py --results /path/to/results
   ```

---

## 🛡 **Benchmarked Tools**

This repository currently supports benchmarking the following tools:

- **Casanovo:** [https://github.com/Noble-Lab/casanovo]
- **pi-HelixNovo:** [https://github.com/PHOENIXcenter/pi-HelixNovo]
- **ContraNovo:** [https://github.com/BEAM-Labs/ContraNovo]

---

## 🤝 **Contributing**

We welcome contributions! To contribute:

1. Fork the repository.
2. Create a new branch for your feature/bugfix.
3. Submit a pull request with a detailed description of your changes.

---


## 📚 **Citation**

If you use this repository in your work, please cite it as follows:

```
@misc{antibody_ms_benchmark,
  author = {Wenbin Jiang, Ling Luo},
  title = {Antibody Mass Spectrometry De Novo Benchmark},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished = {\url{https://github.com/dumbgoos/AbNovoBench}},
}
```

---

## 📜 **License**

This project is licensed under the [MIT License](LICENSE). Feel free to use, modify, and distribute this repository as per the license terms.


