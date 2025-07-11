# Cyanobacterial resource allocation strategies explained by Growth Balance Analysis

This repository contains **R** scripts for analyzing Growth Balance Analysis (GBA) models of *Synechocystis sp. PCC 6803* and generating publication-quality plots. The workflows cover:

- **Base model**
- **Proteome allocation**
- **Extended model**

---

## 📂 Repository Structure

```
Synechocystis_plots/
├── data/                 # (optional) raw and processed data files
├── scripts/              # all R scripts for model analysis and plotting
│   ├── plot_base_model.R
│   ├── plot_proteome.R
│   └── plot_extended_model.R
├── figures/              # output figures (PDF/PNG)
├── README.md             # this file
└── .gitignore            # ignore R history, data caches, etc.
```

---

## 🛠️ Requirements

- **R** (≥ 4.0.0)
- R packages:
  - **ggplot2**
  - **dplyr**
  - **tibble**
  - **writexl** or **openxlsx** (if exporting Excel)

Install missing packages in R:

```r
install.packages(c("ggplot2", "dplyr", "tibble", "writexl", "openxlsx"))
```

---

## ⚙️ Usage

1. **Clone** the repo:
   ```bash
   ```

git clone [https://github.com/Sijr73/Synechocystis\_plots.git](https://github.com/Sijr73/Synechocystis_plots.git) cd Synechocystis\_plots

````

2. **Run** the scripts in R or via command line. For example:
```bash
Rscript GBA.R
afterwards, check the output in Results
````

3. **Customize**: edit any script to point at your own data paths or tweak axis limits, colors, and annotations.

---

## 📈 Scripts Overview

### `plot_base_model.R`

- Reads `newSJM2`, `newSJM2low`, `faizi`, and `faizi_low` model results.
- Plots growth rate vs. % light intensity with error bars for observed data.
- Saves `FigGrowthBase.pdf`.

### `plot_proteome.R`

- Uses `syn_Ribosomefinal` observational data.
- Plots ribosome proteome mass fraction vs. growth rate.
- Saves `FigRibosomeBase.pdf`.

### `plot_extended_model.R`

- Loads `Extendedmodel04032024` and `Extendedmodellow04032024`.
- Plots extended model (with inhibition) and compares to base model and observations.
- Saves `FigGrowthExtended.pdf`.

---

## 🤝 Contributing

1. Fork the repository.
2. Create a new branch: `git checkout -b feature/your-feature`
3. Commit your changes: `git commit -m "Add new analysis script"`
4. Push to your branch: `git push origin feature/your-feature`
5. Open a Pull Request and describe your changes.

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

*Happy plotting!*

