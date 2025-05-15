# Mixed DAIseg

Mixed DAIseg is a Python library for genomic segmentation analysis, implementing algorithms from [DAIseg](https://github.com/Genomics-HSE/DAIseg), [DAIseg.mex](https://github.com/Genomics-HSE/DAIseg.mex), and [DAIseg-general](https://github.com/LeoPlanche/DAIseg). It provides a flexible interface for running segmentation tasks with configurable parameters. The library requires a virtual environment with `numpy` for operation.

---

## Installation

1. **Create a virtual environment** (recommended):
   ```bash
   python3 -m venv daiseg_env
   source daiseg_env/bin/activate  # Linux/MacOS
   # or
   daiseg_env\Scripts\activate.bat  # Windows
   ```

2. **Install dependencies**:
   ```bash
   pip install numpy
   ```

---

## Usage

To run Mixed DAIseg, execute the following command from your working directory:
```bash
python3 daiseg.py --mode <MODE> [FLAGS]
```

### Modes and Required Flags

#### `--mode simple`
Corresponds to the implementation from [DAIseg](https://github.com/Genomics-HSE/DAIseg).  
**Mandatory flags**:

| Flag               | Description                                                                 | Required Value         |
|--------------------|-----------------------------------------------------------------------------|------------------------|
| `--bed`            | Region bed file                                                             | Yes                    |
| `--EM`             | Whether or not to use EM algorithm                                          | Yes                    |
| `--EM_steps`       | Number of EM iterations.                                                    | No                     |
| `--EM_samples`     | Number of samples used in the EM algorithm.                                 | No                     |
| `--HMM_par`        | Path to HMM parameter file.                                                 | Yes                    |
| `--out_prefix`     | Prefix for output files.                                                    | Yes                    |
| `--EM_est`         | Make estimation of the all parameters or only coalescent times              | No                     |
| `--prepared_file`  | Path to a preprocessed file                                                 | Yes                    |
| `--arch_cover`     |                                                                             | Yes                    |
| `--obs_samples`    | File with samples names                                                     | Yes                    |
| `--decoding`       | Viterbi or aposteriory decoding'                                            | Yes                    |
| `--cut_off`        | Decoding cut off                                                            | Yes                    |

#### `--mode admixed`
Corresponds to the implementation from [DAIseg.mex](https://github.com/Genomics-HSE/DAIseg.mex).  
**Mandatory flags**:

| Flag                  | Description                                                                 | Required Value         |
|-----------------------|-----------------------------------------------------------------------------|------------------------|
| `--bed`               | Region bed file                                                             | Yes                    |
| `--EM`                | Whether or not to use EM algorithm                                          | Yes                    |
| `--EM_steps`          | Number of EM iterations.                                                    | No                     |
| `--EM_samples`        | Number of samples used in the EM algorithm.                                 | No                     |
| `--HMM_par`           | Path to HMM parameter file.                                                 | Yes                    |
| `--out_prefix`        | Prefix for output files.                                                    | Yes                    |
| `--prepared_file`     | Path to a preprocessed file                                                 | Yes                    |
| `--arch_cover`        |                                                                             | Yes                    |
| `--obs_samples`       | File with samples names                                                     | Yes                    |
| `--obs_type`          | Type of observed data (e.g., `reads`, `coverage`).                          | No                     |
| `--transition_matrix` | Path to the transition matrix file for HMM.                                 | Yes                    |

#### `--mode general`
Corresponds to the implementation from [DAIseg-general](https://github.com/LeoPlanche/DAIseg).  
**Mandatory flags**:

| Flag            | Description                                                                                 | Required Value         |
|-----------------|---------------------------------------------------------------------------------------------|------------------------|
| `--ind`         | Ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2                     | Yes                    |
| `--vcfOut`      | Path to list of comma-separated vcf/bcf file(s) containing the individuals in the outgroups | Yes                    |
| `--vcfIn`       | Path to list of comma-separated vcf/bcf file(s) containing the individuals in the ingroup   | Yes                    |
| `--dem`         | Demographic model configuration file.                                                       | No                     |
| `--weights`     | file with callability (defaults to all positions being called)                              | No                     |
| `--out_prefix`  | Prefix for output files.                                                                    | No                     |
| `--ancestral`   | fasta file with ancestral information - comma-separated list or wildcards like vcf argument | No                     |
| `--refgenome`   | fasta file with reference genome - comma-separated list or wildcards like vcf argument      | No                     |
| `--haploid`     | Change from using diploid data to haploid data                                              | No                     |
| `--obs_type`    | Use conditional probability                                                                 | No                     |

---

For issues or contributions, refer to the original repositories:  
- [DAIseg](https://github.com/Genomics-HSE/DAIseg)  
- [DAIseg.mex](https://github.com/Genomics-HSE/DAIseg.mex)  
- [DAIseg-general](https://github.com/LeoPlanche/DAIseg)  
```
