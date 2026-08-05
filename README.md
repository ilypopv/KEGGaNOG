# KEGGaNOG

<img src="https://github.com/ilypopv/KEGGaNOG/blob/main/imgs/KaN_logo_light.png#gh-light-mode-only" align="left" width = 25%/>
<img src="https://github.com/ilypopv/KEGGaNOG/blob/main/imgs/KaN_logo_dark.png#gh-dark-mode-only" align="left" width = 25%/>

<br>
<br>

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![KEGG-Decoder](https://img.shields.io/badge/Dependecy-KEGG_Decoder-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Downloads](https://static.pepy.tech/badge/kegganog)](https://pepy.tech/project/kegganog)
[![codecov](https://codecov.io/gh/ilypopv/KEGGaNOG/graph/badge.svg)](https://codecov.io/gh/ilypopv/KEGGaNOG)

![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)
![macOS](https://img.shields.io/badge/mac%20os-000000?style=for-the-badge&logo=macos&logoColor=F0F0F0)

<br>
<br>

## Motivation

[**`eggNOG-mapper`**](https://github.com/eggnogdb/eggnog-mapper) 🤝 [**`KEGG-Decoder`**](https://github.com/bjtully/BioData/blob/master/KEGGDecoder/README.md)

- `eggNOG-mapper` is a comprehensive tool for fast functional annotation of novel sequences. Yet it does not provide any visualization functions.
- `KEGG-Decoder` is a perfect tool for visualizing KEGG Pathways. But it only takes `KEGG-Koala` outputs as an input (including blastKOALA, ghostKOALA, KOFAMSCAN).
- `KEGG-Koala` is a web-tool which can work for more than 24 hours. `eggNOG-mapper` can be installed locally on your PC / server and work faster.
- This tool `KEGGaNOG` makes `eggNOG-mapper` meet `KEGG-Decoder`! It parses `eggNOG-mapper` output, make it fit for the input to `KEGG-Decoder` and then visualize KEGG Pathways as the heatmap!
- **Pro-tip:** `eggNOG-mapper` and `KEGGaNOG` could be wrapped into 🐍 `Snakemake` pipeline making metabolic profiling a "one-click" process!

## Installation

```bash
# Linux / WSL / macOS
conda create -n kegganog pip -y
conda activate kegganog
pip install kegganog
```

## Usage Guide

### Command-line mode

```text
Usage: KEGGaNOG [OPTIONS]                                                      
                                                                                
 KEGGaNOG: Link eggNOG-mapper and KEGG-Decoder for pathway visualization.       
                                                                                
╭─ Options ──────────────────────────────────────────────────────────────────────────────╮
│ --input      -i         TEXT                    Path to eggNOG-mapper annotation file. │
│ --output     -o         TEXT                    Output folder to saveresults.          │
│ --multi      -M                                 Run KEGGaNOG in multi-sample cohort    │
│                                                 profile mode.                          │
│ --overwrite  -overwrite                         Overwrite the output directory if it   │
│                                                 already exists.                        │
│ --dpi        -dpi       INTEGER                 DPI resolution mapping index for the   │
│                                                 output image visualization.            │
│                                                 [default: 300]                         │
│ --color      -c         [Blues|Greens|Reds|     Target seaborn color map palette       │
│                         Purples|Greys|Oranges]  matrix rule.                           │
│                                                 [default: Blues]                       │
│ --name       -n         TEXT                    Sample identity text string for        │
│                                                 axis labeling.                         │
│                                                 [default: SAMPLE]                      │
│ --group      -g                                 Group the pathway matrix heatmap rows  │
│                                                 based on predefined functional         │
│                                                 categories.                            │
│ --web                                           Launch the interactive local web UI    │
│                                                 dashboard at http://localhost:8000.    │
│ --version    -V                                 Show version and exit.                 │
│ --help       -h                                 Show this message and exit.            │
╰────────────────────────────────────────────────────────────────────────────────────────╯
```

🔗 Please visit [KEGGaNOG wiki](https://github.com/ilypopv/KEGGaNOG/wiki) page

### Web interface mode

For an interactive, browser-based experience with live preview and advanced visualization options:

```bash
KEGGaNOG --web
```

Then open `http://localhost:8000` in your browser.

**Features:**

- **Live preview** — visualize plots in real-time as you adjust parameters
- **Interactive settings** — no command-line arguments needed; drop files, tweak colors and dimensions through an intuitive UI
- **Multi-sample analysis** — compare samples using heatmaps, radarplots, correlation networks, stacked bars, and streamgraphs
- **Re-render on the fly** — modify plot parameters without re-running the full analysis (multi mode only)
- **Download results** — export individual plots or the complete results ZIP

## Output examples gallery

### Default visualization

|Single mode|Multi mode|
|-----------|----------|
|![heatmap_figure](https://github.com/user-attachments/assets/2b50518d-1fff-46d6-8bfc-6a5b8c31356d)|![heatmap_figure](https://github.com/user-attachments/assets/484077b9-8212-4aa2-8a3e-1a831afba26f)|

These figures are generated using functional groupping mode (`-g`/`--group`) and `Greens` colormap

### User APIs visualization

|[Barplot](https://github.com/ilypopv/KEGGaNOG/wiki/Barplot-API)|[Boxplot](https://github.com/ilypopv/KEGGaNOG/wiki/Boxplot-API)|[Radarplot](https://github.com/ilypopv/KEGGaNOG/wiki/Radarplot-API)|[Correlation Network](https://github.com/ilypopv/KEGGaNOG/wiki/Correlation-Network-API)|
|-------|-------|---------|-------------------|
|![image](https://github.com/user-attachments/assets/81d69bef-f69c-4960-b2d3-73e348e3853a)|![image](https://github.com/user-attachments/assets/f98fd993-20b7-4b00-b203-83b40fe35f9c)|![image](https://github.com/user-attachments/assets/dd75e5d8-e3c8-4eaa-b009-02c042534a53)|![image](https://github.com/user-attachments/assets/e76057b9-bcfd-4ba9-a4cf-cb7b4269441a)|

|[Stacked Barplot](https://github.com/ilypopv/KEGGaNOG/wiki/Stacked-Barplot-API)|[Streamgraph](https://github.com/ilypopv/KEGGaNOG/wiki/Streamgraph-API)|[Stacked Barplot + Streamgraph](https://github.com/ilypopv/KEGGaNOG/wiki/Combined-Stacked-Barplot-&-Streamgraph)|
|-------|-------|-------|
|![kgnstbar_OLD](https://github.com/user-attachments/assets/11e9e265-52c7-41b7-a284-64f3181caac3)|![kgnstream_OLD](https://github.com/user-attachments/assets/e82654fc-478a-4233-8478-f2c69ee4a1a6)|![combined_white_OLD](https://github.com/user-attachments/assets/6059da3c-4b74-47a2-af1c-427180f44845)|

## Advantages

1. **Seemless Access to KEGG Annotations:** Provides KEGG Ortholog (KO) annotations without requiring a KEGG license.
2. **High-Throughput Capability:** Optimized for rapid KO assignment in large-scale datasets, ideal for metagenomics and genomics projects.
3. **Broad Functional Coverage:** Leverages the extensive eggNOG database to annotate genes across a wide range of taxa.

## Limitation

1. **Indirect KO Mapping:** `eggNOG-mapper` doesn’t directly use the KEGG database, its KO term assignments are inferred through orthologous groups (eggNOG entries). This can sometimes result in less precise annotations.

## Tool name background

`KEGGaNOG` stands for “KEGG out of NOG”, highlighting its purpose: extracting KEGG Ortholog annotations from eggNOG’s Non-supervised Orthologous Groups.

## Citation

If you use `KEGGaNOG` in your research, please cite:

[Popov, I.V., Chikindas, M.L., Venema, K., Ermakov, A.M. and Popov, I.V., 2025.
KEGGaNOG: A Lightweight Tool for KEGG Module Profiling From Orthology-Based Annotations.
Molecular Nutrition & Food Research, p.e70269.
doi.org/10.1002/mnfr.70269](https://onlinelibrary.wiley.com/share/author/QAEWYFJAQ6TAZGDW6FJZ?target=10.1002/mnfr.70269)

## Acknowledgements

For now `KEGGaNOG` uses [**`KEGG-Decoder`**](https://github.com/bjtully/BioData/blob/master/KEGGDecoder/KEGG_decoder.py) as a main dependecy. I greatly thank [**`KEGG-Decoder`**](https://github.com/bjtully/BioData/blob/master/KEGGDecoder/KEGG_decoder.py)'s developers.
