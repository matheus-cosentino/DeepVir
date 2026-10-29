# 🇬🇧 DeepVir Pipeline Structure (English)

Below is the detailed documentation of the 3-directory repository architecture and execution DAGs of the DeepVir pipeline.

## 📊 Pipeline DAGs

### Pipeline Presentation

```mermaid
graph TD;
    A[Raw Reads fastq.gz] --> B(fastp Quality Control);
    B -->|Filtered Reads| C{Analysis Path};
    C -->|Read-Level Taxonomy| D[Kraken2, Diamond & MEGAN];
    C -->|De Novo Assembly| E[SPAdes / Megahit / Flye / Raven];
    
    E -->|Contigs| F(Diamond Annotation vs NR);
    E -->|Contigs| G(Kraken2 Annotation);
    E -->|Contigs| H(Dark Matter Module);
    
    H -->|ORF Prediction & HMM| I(PalmAnnot2);
    H -->|Optional| J(RVDB Screening);
    I -->|Taxonomy Summarization| K[Final RdRp Hits & Annotations];

    D --> L[Krona & MultiQC Report Dashboards];
    F --> L;
    K --> L;
    F --> K;
    J --> K;
```

### Execution DAG (Slurm / Local Profile)

```mermaid
%%{init: {"themeVariables": {"fontSize": "22px"}}}%%
graph TD;
    subgraph Quality_Check ["Quality Check"]
        2["fastp_paired - sample: test"]
    end
    subgraph Diversity_Diamond ["Diversity Diamond"]
        4["append_lineage"]
        5["split_hits_by_taxid"]
        6["map_accession_to_taxid"]
        7["diamond_blastx_contigs - tool: spades_kauto"]
        21["append_lineage"]
        22["split_hits_by_taxid"]
        23["map_accession_to_taxid - tool: reads"]
        24["diamond_blastx_reads"]
    end
    subgraph Diversity_Kraken ["Diversity Kraken"]
        27["kraken_biom_reads - paired: paired"]
        28["kraken2_reads_paired"]
        31["kraken_biom_merge_all"]
        32["kraken_biom_contig"]
        33["kraken2_contigs - tool: spades_kauto"]
    end
    subgraph Dark_Matter_Search ["Dark Matter Search"]
        9["fev2tsv_single"]
        10["palm_annot"]
        11["cd_hit"]
        12["find_orfs"]
        13["get_nohit_fasta"]
        14["dm_validate"]
        15["rvdb_validate_structural"]
        16["rvdb_summarize"]
        17["rvdb_search"]
        20["rvdb_validate_pol"]
    end
    subgraph Assembly ["Assembly"]
        1["spades - kmer_val: auto"]
    end
    subgraph Reporting ["Reporting"]
        3["report_summarize"]
        25["krona_diamond"]
        29["krona_reads_kraken - read_type: paired"]
        30["kraken_rarefaction_plot"]
        34["krona_kraken2"]
        35["krona_diamond"]
    end
    subgraph Setup_Databases ["Setup & Databases"]
        8["download_prot"]
        18["prepare_rvdb"]
        19["get_rvdb"]
        26["krona_update_taxonomy"]
    end
    subgraph Other ["Other"]
        0["all"]
    end
    1 --> 0
    3 --> 0
    14 --> 0
    15 --> 0
    20 --> 0
    21 --> 0
    25 --> 0
    27 --> 0
    29 --> 0
    30 --> 0
    32 --> 0
    34 --> 0
    4 --> 0
    35 --> 0
    2 --> 1
    1 --> 3
    4 --> 3
    9 --> 3
    5 --> 4
    6 --> 5
    7 --> 6
    8 --> 6
    1 --> 7
    10 --> 9
    11 --> 10
    12 --> 11
    13 --> 12
    1 --> 13
    4 --> 13
    1 --> 14
    11 --> 14
    9 --> 14
    1 --> 15
    11 --> 15
    16 --> 15
    17 --> 16
    18 --> 16
    11 --> 17
    18 --> 17
    19 --> 18
    1 --> 20
    11 --> 20
    16 --> 20
    22 --> 21
    23 --> 22
    24 --> 23
    8 --> 23
    2 --> 24
    23 --> 25
    26 --> 25
    8 --> 26
    28 --> 27
    2 --> 28
    28 --> 29
    31 --> 30
    28 --> 31
    33 --> 32
    1 --> 33
    33 --> 34
    6 --> 35
    26 --> 35
```

---

## 📂 Root Directory (`/`) — Launchers & Entry Points

*   **`DeepVir.py`**: Universal pure-Python orchestrator for the entire pipeline. It parses CLI arguments, launches the Desktop GUI (`--gui`), bridges Windows executions to WSL2 automatically, resolves Conda environments, writes dynamic run overrides (`logs/run_overrides.yaml`), and triggers Snakemake.
*   **`DeepVir_gui.bat`**: 1-click launcher for **Windows** (detects native Python or falls back to WSL2).
*   **`DeepVir_gui.command`**: 1-click launcher for **macOS** (Finder double-click).
*   **`DeepVir_gui.sh`**: 1-click launcher for **Linux / WSL**.
*   **`README.md`**: Main project overview, quick-start instructions, and links to tutorials.

---

## 📂 1. `pipeline/` — Pipeline Backend & GUI Engine

Contains the entire codebase, configuration files, Snakemake workflow, and graphical interface:

*   **`pipeline/config/config.yaml`**: Static global workflow parameters (SPAdes K-mer sizes, E-value cutoffs, Flye/Medaka settings, default database paths).
*   **`pipeline/gui/deepvir_gui.py`**: Interactive Tkinter/TTK Desktop Dashboard for configuring parameters, monitoring real-time execution, and opening HTML reports.
*   **`pipeline/profiles/`**: Snakemake computational profiles:
    *   `local/`: For local workstations, laptops, and WSL2 environments.
    *   `profile_slurm/`: For High-Performance Computing (HPC) clusters managed by Slurm.
*   **`pipeline/resources/`**: External databases (`diamond/`, `kraken2/`, `taxonomy/`), PalmAnnot HMM models (`palm_annot/`), and institutional logos (`logo/`).
*   **`pipeline/workflow/`**:
    *   `Snakefile`: Master Snakemake entry point.
    *   `multiqc_config.yaml`: Visual layout configuration for MultiQC reports.
    *   `envs/`: Reproducible Conda environment YAML specifications for each bioinformatics tool (`spades`, `diamond`, `kraken2`, `kraken2_biom`, `palm_annot`, `report`, etc.).
    *   `rules/`: Modular Snakemake rule files (`.smk`) including `get_data.smk`, `fastp.smk`, `denovo.smk`, `diamond.smk`, `kraken2.smk`, `lineage.smk`, `megan.smk`, `basta.smk`, `darkmatter.smk`, `RVDB.smk`, `reporting.smk`, and `krona.smk`.
    *   `scripts/`: Helper Python and R scripts (`functions.py`, `dm_validate.py`, `palm_annot_run.py`, `generate_report.R`, `Report_Model.Rmd`).

---

## 📂 2. `docs/` — License & Documentation

Contains the software license and all Markdown tutorials:
*   **`LICENSE`**: GNU General Public License v3.0 (GPL v3).
*   **`local_use.md`**: First-time local setup and Python CLI / GUI quick start.
*   **`desktop_gui.md`**: Complete guide to the Desktop GUI and visual components.
*   **`full_databases_setup.md`**: Production database preparation for Kraken2 and DIAMOND NR.
*   **`hpc_use.md`**: Running DeepVir on Slurm HPC clusters with `python3 DeepVir.py`.
*   **`advanced_tutorial.md`**: Tuning `fastp`, assemblers, and downstream parameters.
*   **`public_data_usage.md`**: Downloading and analyzing NCBI SRA datasets.
*   **`expected_outputs.md`**: Overview of `results/` and `logs/` structures.
*   **`structure_and_dags/workflow_structure.md`**: This trilingual architecture and DAG guide.

---

## 📂 3. `logs/` — Centralized Execution Logs

Stores all runtime logs and configuration snapshots:
*   **`logs/<sample>/*.log`**: Per-sample logs generated by each Snakemake rule.
*   **`logs/workflow/*.log`**: Global workflow logs (MultiQC, Krona taxonomy, RVDB database setup).
*   **`logs/run_overrides.yaml`**: Dynamic configuration generated by `DeepVir.py` for the active run (plus timestamped archives `run_overrides_<timestamp>.yaml`).
*   **`logs/deepvir_desktop_gui.log`**: Diagnostic log from the 1-click GUI launchers.

<br>
<hr>
<br>

# 🇧🇷 Estrutura do Pipeline DeepVir (Português)

Abaixo está a documentação detalhada da arquitetura em 3 diretórios principais e dos pontos de entrada em Python do pipeline DeepVir.

## 📂 Diretório Raiz (`/`) — Executores e Atalhos 1-Clique

*   **`DeepVir.py`**: Executor universal em Python puro do pipeline. Recebe parâmetros via linha de comando (CLI) ou abre a interface gráfica (`--gui`), traduz caminhos automaticamente entre Windows e WSL2, gerencia os ambientes Conda, gera `logs/run_overrides.yaml` e aciona o Snakemake.
*   **`DeepVir_gui.bat`**: Executável de 1 clique para **Windows** (usa Python nativo ou WSL2 automaticamente).
*   **`DeepVir_gui.command`**: Executável de 1 clique para **macOS** (duplo-clique no Finder).
*   **`DeepVir_gui.sh`**: Executável de 1 clique para **Linux / WSL**.
*   **`README.md`**: Documentação principal do projeto, instruções rápidas de uso em Python e links para os tutoriais.

---

## 📂 1. `pipeline/` — Código-Fonte do Pipeline e Backend

Concentra todo o motor do fluxo de trabalho, configurações e interface gráfica:

*   **`pipeline/config/config.yaml`**: Parâmetros técnicos globais (K-mers do SPAdes, limiares de E-value, modelos Flye/Medaka e caminhos padrão).
*   **`pipeline/gui/deepvir_gui.py`**: Painel gráfico Desktop (Tkinter/TTK) para configuração, monitoramento em tempo real e visualização de relatórios HTML.
*   **`pipeline/profiles/`**: Perfis computacionais do Snakemake (`local/` para execução local/WSL e `profile_slurm/` para clusters HPC com Slurm).
*   **`pipeline/resources/`**: Bancos de dados (`diamond/`, `kraken2/`, `taxonomy/`), modelos HMM do PalmAnnot (`palm_annot/`) e logotipos (`logo/`).
*   **`pipeline/workflow/`**:
    *   `Snakefile`: Arquivo mestre de entrada do Snakemake.
    *   `multiqc_config.yaml`: Configurações visuais do relatório MultiQC.
    *   `envs/`: Arquivos YAML dos ambientes Conda reprodutíveis para cada ferramenta.
    *   `rules/`: Regras modulares do Snakemake (`.smk`: `get_data.smk`, `fastp.smk`, `denovo.smk`, `diamond.smk`, `kraken2.smk`, `lineage.smk`, `megan.smk`, `basta.smk`, `darkmatter.smk`, `RVDB.smk`, `reporting.smk`, `krona.smk`).
    *   `scripts/`: Scripts auxiliares em Python e R (`functions.py`, `dm_validate.py`, `palm_annot_run.py`, `generate_report.R`, `Report_Model.Rmd`).

---

## 📂 2. `docs/` — Licenças e Documentações

Reúne a licença legal e todos os manuais em Markdown:
*   **`LICENSE`**: Termos legais da licença GNU GPL v3.
*   **`local_use.md`**, **`desktop_gui.md`**, **`full_databases_setup.md`**, **`hpc_use.md`**, **`advanced_tutorial.md`**, **`public_data_usage.md`**, **`expected_outputs.md`** e **`structure_and_dags/workflow_structure.md`**.

---

## 📂 3. `logs/` — Logs Centralizados de Execução

Armazena todos os registros de execução e snapshots de configuração:
*   **`logs/<sample>/*.log`**: Logs individuais de cada regra do Snakemake por amostra.
*   **`logs/workflow/*.log`**: Logs globais do workflow (MultiQC, taxonomia Krona, preparação RVDB).
*   **`logs/run_overrides.yaml`**: Arquivo dinâmico gerado pelo `DeepVir.py` a cada rodada (com histórico `run_overrides_<timestamp>.yaml`).
*   **`logs/deepvir_desktop_gui.log`**: Log de inicialização da interface gráfica Desktop.

<br>
<hr>
<br>

# 🇫🇷 Structure du Pipeline DeepVir (Français)

Vous trouverez ci-dessous la documentation détaillée de l'architecture en 3 répertoires principaux du pipeline DeepVir et de son exécuteur universel en Python.

## 📂 Répertoire Racine (`/`) — Exécuteur Python & Lanceurs 1-Clic

*   **`DeepVir.py`**: Orchestrateur universel en Python pur. Il gère les arguments CLI, lance l'interface graphique (`--gui`), assure le pont automatique Windows ↔ WSL2, active Conda, génère `logs/run_overrides.yaml` et déclenche Snakemake.
*   **`DeepVir_gui.bat`**: Lanceur 1-clic pour **Windows** (avec support WSL2).
*   **`DeepVir_gui.command`**: Lanceur 1-clic pour **macOS** (Finder).
*   **`DeepVir_gui.sh`**: Lanceur 1-clic pour **Linux / WSL**.
*   **`README.md`**: Documentation principale du projet.

---

## 📂 1. `pipeline/` — Code Backend & Interface Graphique

*   **`pipeline/config/config.yaml`**: Paramètres techniques globaux (tailles de K-mers SPAdes, seuils E-value, options Flye/Medaka).
*   **`pipeline/gui/deepvir_gui.py`**: Tableau de bord graphique Desktop (Tkinter/TTK).
*   **`pipeline/profiles/`**: Profils d'exécution Snakemake (`local/` et `profile_slurm/`).
*   **`pipeline/resources/`**: Bases de données, modèles HMM (`palm_annot/`) et logos (`logo/`).
*   **`pipeline/workflow/`**: Contient le `Snakefile`, `multiqc_config.yaml`, les environnements Conda (`envs/`), les règles modulaires (`rules/*.smk`) et les scripts Python/R (`scripts/`).

---

## 📂 2. `docs/` — Licence & Documentation Markdown

*   **`LICENSE`**: Licence logicielle GNU GPL v3.
*   Tutoriels complets en Markdown (`local_use.md`, `desktop_gui.md`, `full_databases_setup.md`, `hpc_use.md`, `advanced_tutorial.md`, `public_data_usage.md`, `expected_outputs.md`, `structure_and_dags/workflow_structure.md`).

---

## 📂 3. `logs/` — Journaux d'Exécution Centralisés

*   **`logs/<sample>/*.log`**: Journaux d'exécution Snakemake par échantillon et par règle.
*   **`logs/workflow/*.log`**: Journaux globaux du workflow.
*   **`logs/run_overrides.yaml`**: Configuration dynamique générée par `DeepVir.py` à chaque exécution.
*   **`logs/deepvir_desktop_gui.log`**: Journal de diagnostic de l'interface graphique.
