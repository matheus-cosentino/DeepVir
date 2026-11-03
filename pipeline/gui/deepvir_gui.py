#!/usr/bin/env python3
"""
DeepVir GUI - Desktop Graphical User Interface
=============================================================================
Author: MSc. Matheus Cosentino / DeepVir Team
Description: Modern Desktop GUI in Python for configuring, running, and monitoring
the DeepVir viral metagenomics pipeline both Locally (WSL/Linux/macOS) and on
Remote High-Performance Computing Clusters (Slurm/HPC).
=============================================================================
"""

import os
import re
import sys
import json
import gzip
import base64
import subprocess
import threading
import queue
import time
import webbrowser
import platform
import html
from html.parser import HTMLParser
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, font as tkfont

# Regular expression to match ANSI escape sequences (colors and cursor control)
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]|\r")

# =============================================================================
# CONSTANTS & GLOBAL PATH CONFIGURATION
# =============================================================================
APP_TITLE = "DeepVir 2.0 - Viral Metagenomics & Dark Matter Discovery"
APP_VERSION = "v2.0 Desktop"

# Detect 3-directory project layout relative to this script (pipeline/gui/deepvir_gui.py)
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
PROJECT_DIR = PIPELINE_DIR.parent if (PIPELINE_DIR.parent / "DeepVir.py").exists() else PIPELINE_DIR
DOCS_DIR = PROJECT_DIR / "docs"
LOGS_DIR = PROJECT_DIR / "logs"
PROFILES_DIR = PIPELINE_DIR / "profiles"

# Color Theme: Modern Dark Slate & Emerald Aesthetic
COLOR_BG = "#0f172a"            # Slate 900 (Main Window Background)
COLOR_CARD = "#1e293b"          # Slate 800 (Card Panels & Containers)
COLOR_CARD_BORDER = "#334155"   # Slate 700 (Subtle Borders)
COLOR_TEXT = "#f8fafc"          # Slate 50 (High-contrast Primary Text)
COLOR_MUTED = "#cbd5e1"         # Slate 300 (Readable Secondary / Help Text)
COLOR_PRIMARY = "#10b981"       # Emerald 500 (Primary Action Buttons & Highlights)
COLOR_PRIMARY_HOVER = "#059669" # Emerald 600 (Button Hover State)
COLOR_ACCENT = "#38bdf8"        # Sky 400 (Section Headers & Accents)
COLOR_DANGER = "#ef4444"        # Red 500 (Stop / Abort Button)
COLOR_CONSOLE_BG = "#0b0f19"    # Deep Dark (Terminal Background)
COLOR_CONSOLE_TEXT = "#34d399"  # Light Emerald (Terminal Monospace Font)

# Discrete zoom presets compatible with BOTH X11 bitmap fonts (8, 10, 12, 14, 18, 24)
# AND TrueType / FreeType / DirectWrite vector fonts on Windows, macOS, and Linux.
ZOOM_PRESETS = [
    {"label": "80%",  "small": 8,  "body": 10, "subhead": 12, "head": 12, "title": 14, "mono": 10},
    {"label": "90%",  "small": 10, "body": 10, "subhead": 12, "head": 14, "title": 18, "mono": 10},
    {"label": "100%", "small": 10, "body": 12, "subhead": 14, "head": 14, "title": 18, "mono": 12},
    {"label": "125%", "small": 12, "body": 14, "subhead": 18, "head": 18, "title": 24, "mono": 14},
    {"label": "150%", "small": 14, "body": 18, "subhead": 18, "head": 24, "title": 24, "mono": 18},
    {"label": "180%", "small": 18, "body": 24, "subhead": 24, "head": 24, "title": 24, "mono": 24},
]
DEFAULT_ZOOM_INDEX = 2  # Starts at 100% (12pt body / 14pt headers / 18pt title)

# Canonical rules and their descriptions for the Per-Rule Resource Editor
RULE_DEFINITIONS = [
    ("bwa_index_host", "BWA Index Host Genome", "Host Filtering"),
    ("host_filter_paired", "Host Read Depletion via BWA (Paired)", "Host Filtering"),
    ("host_filter_unpaired", "Host Read Depletion via BWA (Unpaired)", "Host Filtering"),
    ("spades", "De Novo Assembly (SPAdes)", "Assembly"),
    ("megahit", "De Novo Assembly (MEGAHIT)", "Assembly"),
    ("flye", "De Novo Assembly (Flye Long-Read)", "Assembly"),
    ("diamond_blastx_contigs", "DIAMOND BLASTx on Assembled Contigs", "DIAMOND & Taxonomy"),
    ("diamond_blastx_reads", "DIAMOND BLASTx on Filtered Reads", "DIAMOND & Taxonomy"),
    ("reads_diamond_daa", "DIAMOND DAA Alignment for MEGAN", "DIAMOND & Taxonomy"),
    ("map_accession_to_taxid", "NCBI Accession-to-TaxID Mapping", "DIAMOND & Taxonomy"),
    ("basta_search", "BASTA LCA Taxonomic Classification", "DIAMOND & Taxonomy"),
    ("kraken2_contigs", "Kraken2 Classification on Contigs", "Kraken2"),
    ("kraken2_reads_paired", "Kraken2 Classification on Paired Reads", "Kraken2"),
    ("kraken2_reads_unpaired", "Kraken2 Classification on Unpaired Reads", "Kraken2"),
    ("cd_hit", "CD-HIT ORF Clustering (Dark Matter)", "Dark Matter & Reports"),
    ("rvdb_search", "HMMER Search vs RVDB Viral Profiles", "Dark Matter & Reports"),
    ("krona_update_taxonomy", "Krona Taxonomy Database Update", "Dark Matter & Reports"),
]


def _family_actually_scales(root: tk.Tk, family: str) -> bool:
    """
    Tests whether `family` genuinely changes rendered pixel height between 10pt and 18pt
    on the current display server.
    """
    if not family:
        return False
    try:
        f_small = tkfont.Font(root=root, family=family, size=10)
        f_large = tkfont.Font(root=root, family=family, size=18)
        return (
            f_large.metrics("linespace") > f_small.metrics("linespace")
            or f_large.measure("DeepVir") > f_small.measure("DeepVir")
        )
    except Exception:
        return False


def select_ui_font_family(root: tk.Tk) -> str:
    """
    Selects a verified scalable UI font family on Windows, macOS, or Linux/WSL X11.
    """
    try:
        available_list = list(tkfont.families(root))
        available = {f.lower(): f for f in available_list}
    except Exception:
        available_list = []
        available = {}

    try:
        default_actual = tkfont.nametofont("TkDefaultFont").actual("family")
    except Exception:
        default_actual = "helvetica"

    preferred = [
        "Segoe UI",
        "DejaVu Sans",
        "Ubuntu",
        "Liberation Sans",
        "Noto Sans",
        "Helvetica Neue",
        "Arial",
        "Helvetica",
        "helvetica",
        "lucida",
        "times",
        "courier",
        default_actual,
    ]

    for cand in preferred:
        actual_name = available.get(cand.lower(), cand)
        if _family_actually_scales(root, actual_name):
            return actual_name

    for fam in available_list:
        if _family_actually_scales(root, fam):
            return fam

    return default_actual or "helvetica"


def select_mono_font_family(root: tk.Tk) -> str:
    """
    Selects a verified scalable monospace font family for the terminal console.
    """
    try:
        available_list = list(tkfont.families(root))
        available = {f.lower(): f for f in available_list}
    except Exception:
        available_list = []
        available = {}

    try:
        fixed_actual = tkfont.nametofont("TkFixedFont").actual("family")
    except Exception:
        fixed_actual = "courier"

    preferred = [
        "Consolas",
        "DejaVu Sans Mono",
        "Ubuntu Mono",
        "Liberation Mono",
        "Menlo",
        "Courier New",
        "Courier",
        "courier",
        "lucida",
        fixed_actual,
    ]

    for cand in preferred:
        actual_name = available.get(cand.lower(), cand)
        if _family_actually_scales(root, actual_name):
            return actual_name

    for fam in available_list:
        if _family_actually_scales(root, fam):
            return fam

    return fixed_actual or "courier"


def parse_profile_yaml(profile_path: Path) -> dict:
    """
    Lightweight parser for Snakemake profile `config.yaml` (`local` or `profile_slurm`).
    Extracts top-level settings, `default-resources`, and `set-resources` per rule
    without requiring PyYAML in the GUI Python environment.
    """
    data = {
        "jobs": 8,
        "restart_times": 2,
        "latency_wait": 600,
        "runtime": "1-00:00:00",
        "is_slurm": False,
        "default_threads": 1,
        "default_mem_mb": 2000,
        "rules": {},
    }
    for rule_name, _, _ in RULE_DEFINITIONS:
        data["rules"][rule_name] = {"threads": 6, "mem_mb": 12000}

    if not profile_path.is_file():
        return data

    text = profile_path.read_text(encoding="utf-8", errors="replace")
    if "executor: slurm" in text or "profile_slurm" in str(profile_path):
        data["is_slurm"] = True
        data["default_mem_mb"] = 4000
        for rule_name, _, _ in RULE_DEFINITIONS:
            data["rules"][rule_name] = {"threads": 12, "mem_mb": 60000}

    section = None
    current_rule = None

    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()

        if indent == 0:
            current_rule = None
            if stripped.startswith("default-resources:"):
                section = "default-resources"
                continue
            elif stripped.startswith("set-resources:"):
                section = "set-resources"
                continue
            else:
                section = None
                if ":" in stripped:
                    k, v = [p.strip() for p in stripped.split(":", 1)]
                    if k == "jobs" and v.isdigit():
                        data["jobs"] = int(v)
                    elif k == "restart-times" and v.isdigit():
                        data["restart_times"] = int(v)
                    elif k == "latency-wait" and v.isdigit():
                        data["latency_wait"] = int(v)
                    elif k == "runtime" and v:
                        data["runtime"] = v.strip("'\"")
                    elif k == "executor" and "slurm" in v:
                        data["is_slurm"] = True

        elif section == "default-resources" and indent >= 2:
            if ":" in stripped:
                k, v = [p.strip() for p in stripped.split(":", 1)]
                if k == "threads" and v.isdigit():
                    data["default_threads"] = int(v)
                elif k == "mem_mb" and v.isdigit():
                    data["default_mem_mb"] = int(v)

        elif section == "set-resources":
            if indent == 2 and stripped.endswith(":"):
                current_rule = stripped[:-1].strip()
                if current_rule not in data["rules"]:
                    data["rules"][current_rule] = {"threads": 6, "mem_mb": 12000}
            elif indent >= 4 and current_rule and ":" in stripped:
                k, v = [p.strip() for p in stripped.split(":", 1)]
                if k in ("threads", "cpus_per_task") and v.isdigit():
                    data["rules"][current_rule]["threads"] = int(v)
                elif k == "mem_mb" and v.isdigit():
                    data["rules"][current_rule]["mem_mb"] = int(v)

    return data


def write_profile_yaml(profile_path: Path, data: dict) -> None:
    """
    Serializes the profile settings (`default-resources` and `set-resources` per rule)
    back to `pipeline/profiles/<profile>/config.yaml`.
    """
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    is_slurm = data.get("is_slurm", False)
    title_str = "Unified Slurm Profile" if is_slurm else "Local Profile"

    lines = [
        "###################################################################################",
        f"#                              {title_str:<43}#",
        "#                         MSc. Matheus Cosentino                                  #",
        "#                         version: 2.0  |  10.2026                                #",
        "###################################################################################",
        "",
        "######################",
        "# --- 1. Params --- #",
        "#####################",
        f"jobs: {int(data.get('jobs', 8))}",
    ]

    if is_slurm:
        lines.extend([
            "executor: slurm",
            f"restart-times: {int(data.get('restart_times', 0))}",
            f"latency-wait: {int(data.get('latency_wait', 600))}",
            "keep-going: true",
            "rerun-incomplete: true",
            "use-conda: true",
            "show-failed-logs: true",
            f"runtime: {data.get('runtime', '1-00:00:00')}",
        ])
    else:
        lines.extend([
            f"restart-times: {int(data.get('restart_times', 2))}",
            "keep-going: true",
            "rerun-incomplete: true",
            "use-conda: true",
            "show-failed-logs: true",
        ])

    lines.extend([
        "",
        "#######################################",
        "# --- 2. Default Generic Resources --- #",
        "#######################################",
        "default-resources:",
        f"  threads: {int(data.get('default_threads', 1))}",
    ])
    if is_slurm:
        lines.append(f"  cpus_per_task: {int(data.get('default_threads', 1))}")
    lines.append(f"  mem_mb: {int(data.get('default_mem_mb', 2000))}")

    lines.extend([
        "",
        "#######################################",
        "# --- 3. Specific Rule Resources  --- #",
        "#######################################",
        "set-resources:",
    ])

    rules_dict = data.get("rules", {})
    for rule_name, _, _ in RULE_DEFINITIONS:
        r_cfg = rules_dict.get(rule_name, {"threads": 6, "mem_mb": 12000})
        t_val = int(r_cfg.get("threads", 6))
        m_val = int(r_cfg.get("mem_mb", 12000))
        lines.append(f"  {rule_name}:")
        lines.append(f"    mem_mb: {m_val}")
        lines.append(f"    threads: {t_val}")
        if is_slurm:
            lines.append(f"    cpus_per_task: {t_val}")

    profile_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


# =============================================================================
# MAIN APPLICATION CLASS
# =============================================================================
class DeepVirApp(tk.Tk):
    """
    Main DeepVir Desktop Application Window.
    Coordinates all UI components, dynamic font scaling, per-rule resource profiles,
    direct SRA accession entry, threading, and inter-process communication with `DeepVir.py`.
    """

    def __init__(self):
        super().__init__()

        # Configure window properties
        self.title(f"{APP_TITLE} ({APP_VERSION})")
        self.geometry("1320x900")
        self.minsize(1040, 750)
        self.configure(bg=COLOR_BG)

        # Detect verified scalable OS fonts
        self.ui_font = select_ui_font_family(self)
        self.mono_font = select_mono_font_family(self)

        # Current index into ZOOM_PRESETS
        self.zoom_index = DEFAULT_ZOOM_INDEX
        self.zoom_label_var = tk.StringVar(value=ZOOM_PRESETS[self.zoom_index]["label"])

        # Track plain tk widgets for live font rescaling: list of (widget, tier, weight)
        self._scalable_widgets = []
        self._config_canvas = None
        self._profile_canvas = None

        # Per-rule resource editor state variables: {rule_name: (IntVar(threads), IntVar(mem_mb), Label(gb_text))}
        self.rule_resource_vars = {}

        # Process Management State Variables
        self.current_process = None
        self.is_running = False
        self.log_queue = queue.Queue()
        self.active_log_path = None

        # Configure TTK styles and widget themes
        self.setup_styles()

        # Build all visual components, tabs, cards, and buttons
        self.create_widgets()

        # Load initial profile resources from disk (`local`)
        self.load_profile_into_editor("local", silent=True)

        # Apply initial font preset across all widgets
        self.apply_font_preset()

        # Keyboard shortcuts for zooming font size (Ctrl + / Ctrl -)
        self.bind_all("<Control-plus>", lambda e: self.change_font_size(1))
        self.bind_all("<Control-equal>", lambda e: self.change_font_size(1))
        self.bind_all("<Control-minus>", lambda e: self.change_font_size(-1))

        # Start the non-blocking log queue processor (runs every 100 milliseconds)
        self.after(100, self.process_log_queue)

    def _tier_font(self, tier: str, weight: str = "normal"):
        """Returns a `(family, size, weight)` tuple for the given typography tier."""
        preset = ZOOM_PRESETS[self.zoom_index]
        size = preset.get(tier, preset["body"])
        family = self.mono_font if tier == "mono" else self.ui_font
        if weight and weight != "normal":
            return (family, size, weight)
        return (family, size)

    def register_scalable(self, widget, tier: str = "body", weight: str = "normal"):
        """Registers a widget with a typography tier ('small', 'body', 'subhead', 'head', 'title', 'mono')."""
        self._scalable_widgets.append((widget, tier, weight))
        try:
            widget.configure(font=self._tier_font(tier, weight))
        except tk.TclError:
            pass

    def change_font_size(self, step: int):
        """Steps up or down through ZOOM_PRESETS and immediately updates all widgets."""
        new_idx = max(0, min(len(ZOOM_PRESETS) - 1, self.zoom_index + step))
        if new_idx == self.zoom_index:
            return
        self.zoom_index = new_idx
        self.apply_font_preset()

    def apply_font_preset(self):
        """Applies the active ZOOM_PRESETS entry to all named fonts, TTK styles, and tk widgets."""
        preset = ZOOM_PRESETS[self.zoom_index]
        self.zoom_label_var.set(preset["label"])

        # 1. Update Tk's built-in named fonts
        for named, tier, weight in [
            ("TkDefaultFont", "body", "normal"),
            ("TkTextFont", "body", "normal"),
            ("TkMenuFont", "body", "normal"),
            ("TkHeadingFont", "subhead", "bold"),
            ("TkCaptionFont", "body", "bold"),
            ("TkSmallCaptionFont", "small", "normal"),
            ("TkFixedFont", "mono", "normal"),
        ]:
            try:
                nf = tkfont.nametofont(named)
                nf.configure(
                    family=self.mono_font if tier == "mono" else self.ui_font,
                    size=preset[tier],
                    weight=weight,
                )
            except Exception:
                pass

        # 2. Reconfigure TTK styles
        self.setup_styles()

        # 3. Reconfigure all registered tk / ttk widgets
        for widget, tier, weight in self._scalable_widgets:
            try:
                widget.configure(font=self._tier_font(tier, weight))
            except tk.TclError:
                pass

        # 4. Refresh layout and scroll regions immediately
        self.update_idletasks()
        for cv in (self._config_canvas, self._profile_canvas):
            if cv is not None:
                try:
                    cv.configure(scrollregion=cv.bbox("all"))
                except tk.TclError:
                    pass

    def setup_styles(self):
        """
        Configures TTK widget styles using the current zoom preset's font sizes.
        """
        self.style = ttk.Style(self)
        self.style.theme_use("clam")

        # Global base styling
        self.style.configure(".", background=COLOR_BG, foreground=COLOR_TEXT, font=self._tier_font("body"))
        self.style.configure("TFrame", background=COLOR_BG)
        self.style.configure("Card.TFrame", background=COLOR_CARD, relief="flat")

        # Notebook (Tab container) styling
        self.style.configure("TNotebook", background=COLOR_BG, borderwidth=0)
        self.style.configure(
            "TNotebook.Tab",
            background=COLOR_CARD,
            foreground=COLOR_TEXT,
            padding=[18, 10],
            font=self._tier_font("subhead", "bold"),
            borderwidth=0,
        )
        self.style.map(
            "TNotebook.Tab",
            background=[("selected", COLOR_PRIMARY), ("active", "#2d3748")],
            foreground=[("selected", "#ffffff"), ("active", "#ffffff")],
        )

        # Typography and Label variants
        self.style.configure("Header.TLabel", font=self._tier_font("head", "bold"), foreground=COLOR_TEXT, background=COLOR_CARD)
        self.style.configure("Sub.TLabel", font=self._tier_font("small"), foreground=COLOR_MUTED, background=COLOR_CARD)
        self.style.configure("Card.TLabel", background=COLOR_CARD, foreground=COLOR_TEXT, font=self._tier_font("body"))
        self.style.configure("Status.TLabel", font=self._tier_font("body", "bold"), padding=[12, 5])

        # Checkbutton styles
        self.style.configure(
            "TCheckbutton",
            background=COLOR_CARD,
            foreground=COLOR_TEXT,
            font=self._tier_font("body"),
            indicatormargin=[0, 3, 8, 3],
        )
        self.style.map(
            "TCheckbutton",
            background=[("active", COLOR_CARD)],
            foreground=[("active", COLOR_PRIMARY), ("disabled", "#64748b")],
        )

        # Radiobutton styles
        self.style.configure("TRadiobutton", background=COLOR_CARD, foreground=COLOR_TEXT, font=self._tier_font("body", "bold"))
        self.style.map("TRadiobutton", background=[("active", COLOR_CARD)], foreground=[("active", COLOR_ACCENT)])

        # Combobox styling
        self.style.configure("TCombobox", font=self._tier_font("body"), padding=4)
        self.option_add("*TCombobox*Listbox.font", self._tier_font("body"))

        # Primary Action Button (Run Pipeline)
        self.style.configure(
            "Primary.TButton",
            font=self._tier_font("subhead", "bold"),
            background=COLOR_PRIMARY,
            foreground="#ffffff",
            padding=[18, 10],
            borderwidth=0,
        )
        self.style.map("Primary.TButton", background=[("active", COLOR_PRIMARY_HOVER), ("disabled", "#475569")])

        # Danger Button (Stop / Terminate Process)
        self.style.configure(
            "Danger.TButton",
            font=self._tier_font("subhead", "bold"),
            background=COLOR_DANGER,
            foreground="#ffffff",
            padding=[18, 10],
            borderwidth=0,
        )
        self.style.map("Danger.TButton", background=[("active", "#dc2626"), ("disabled", "#475569")])

        # Secondary Button (Browse, Clear, Save, Dry-Run)
        self.style.configure(
            "Secondary.TButton",
            font=self._tier_font("body", "bold"),
            background="#334155",
            foreground=COLOR_TEXT,
            padding=[12, 8],
            borderwidth=0,
        )
        self.style.map("Secondary.TButton", background=[("active", "#475569")])

    def create_widgets(self):
        """
        Constructs the primary UI layout:
        1. Top header banner with title, live font zoom controls (A- / % / A+), mode selector, and status badge.
        2. Main Tabbed Notebook containing:
           - 1. Configuration & Execution
           - 2. Profile & Rule Resources
           - 3. Real-Time Terminal Console
           - 4. Results & Reports
        """
        # ---------------------------------------------------------------------
        # Top Header Banner
        # ---------------------------------------------------------------------
        header_frame = tk.Frame(self, bg=COLOR_CARD, height=72, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        header_frame.pack(fill="x", padx=12, pady=(10, 6))

        # Title branding (100% pure ASCII)
        title_box = tk.Frame(header_frame, bg=COLOR_CARD)
        title_box.pack(side="left", padx=16, pady=12)

        lbl_title = tk.Label(title_box, text="DeepVir 2.0", fg=COLOR_PRIMARY, bg=COLOR_CARD)
        self.register_scalable(lbl_title, "title", "bold")
        lbl_title.pack(side="left")

        lbl_desc = tk.Label(title_box, text=" | Viral Discovery Pipeline Dashboard", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(lbl_desc, "subhead")
        lbl_desc.pack(side="left")

        # Execution Status Badge (pure ASCII)
        self.status_var = tk.StringVar(value="[ READY ]")
        self.lbl_status_badge = tk.Label(
            header_frame,
            textvariable=self.status_var,
            fg="#10b981",
            bg="#064e3b",
            padx=14,
            pady=6,
            relief="flat",
        )
        self.register_scalable(self.lbl_status_badge, "body", "bold")
        self.lbl_status_badge.pack(side="right", padx=16, pady=12)

        # Font Size Zoom Controls (A- / 100% / A+)
        font_box = tk.Frame(header_frame, bg=COLOR_CARD)
        font_box.pack(side="right", padx=10)

        lbl_font = tk.Label(font_box, text="Font:", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(lbl_font, "body")
        lbl_font.pack(side="left", padx=(0, 4))

        btn_font_dec = tk.Button(
            font_box,
            text=" [ A- ] ",
            bg="#334155",
            fg=COLOR_TEXT,
            relief="flat",
            padx=6,
            pady=3,
            cursor="hand2",
            activebackground="#475569",
            activeforeground="#ffffff",
            command=lambda: self.change_font_size(-1),
        )
        self.register_scalable(btn_font_dec, "body", "bold")
        btn_font_dec.pack(side="left", padx=2)

        lbl_zoom_pct = tk.Label(font_box, textvariable=self.zoom_label_var, fg=COLOR_ACCENT, bg=COLOR_CARD, width=5)
        self.register_scalable(lbl_zoom_pct, "body", "bold")
        lbl_zoom_pct.pack(side="left", padx=2)

        btn_font_inc = tk.Button(
            font_box,
            text=" [ A+ ] ",
            bg="#334155",
            fg=COLOR_TEXT,
            relief="flat",
            padx=6,
            pady=3,
            cursor="hand2",
            activebackground="#475569",
            activeforeground="#ffffff",
            command=lambda: self.change_font_size(1),
        )
        self.register_scalable(btn_font_inc, "body", "bold")
        btn_font_inc.pack(side="left", padx=2)

        # Execution Mode Selector (Local PC vs Remote Slurm Cluster)
        mode_box = tk.Frame(header_frame, bg=COLOR_CARD)
        mode_box.pack(side="right", padx=14)

        self.mode_var = tk.StringVar(value="local")
        rb_local = ttk.Radiobutton(
            mode_box,
            text="[Local] WSL / PC",
            value="local",
            variable=self.mode_var,
            command=self.on_mode_change,
        )
        rb_local.pack(side="left", padx=8)

        rb_slurm = ttk.Radiobutton(
            mode_box,
            text="[HPC] Cluster (Slurm)",
            value="slurm",
            variable=self.mode_var,
            command=self.on_mode_change,
        )
        rb_slurm.pack(side="left", padx=8)

        # ---------------------------------------------------------------------
        # Main Tabbed Container
        # ---------------------------------------------------------------------
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=12, pady=6)

        # Tab 1: Configuration & Execution
        self.tab_config = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_config, text="1. Configuration & Execution")
        self.build_tab_config(self.tab_config)

        # Tab 2: Profile & Per-Rule Resources Editor
        self.tab_profile = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_profile, text="2. Profile & Rule Resources")
        self.build_tab_profile(self.tab_profile)

        # Tab 3: Terminal / Live Console
        self.tab_terminal = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_terminal, text="3. Real-Time Terminal Console")
        self.build_tab_terminal(self.tab_terminal)

        # Tab 4: Results & Reports Viewer
        self.tab_results = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_results, text="4. Results & Reports")
        self.build_tab_results(self.tab_results)

        # Automatically refresh the Results file list whenever the user switches to Tab 4
        self.notebook.bind("<<NotebookTabChanged>>", self._on_notebook_tab_changed)

    # =========================================================================
    # TAB 1: CONFIGURATION & EXECUTION BUILDER
    # =========================================================================
    def build_tab_config(self, parent):
        """
        Builds the contents of the Configuration tab with a scrollable content area,
        direct SRA accession input box, and a docked bottom execution bar.
        """
        # Docked Bottom Action Execution Bar (always visible)
        action_bar = tk.Frame(parent, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        action_bar.pack(side="bottom", fill="x", padx=6, pady=(4, 6))

        btn_box = tk.Frame(action_bar, bg=COLOR_CARD)
        btn_box.pack(fill="x", padx=16, pady=12)

        self.btn_dryrun = ttk.Button(
            btn_box,
            text="[Preview] Dry-Run",
            style="Secondary.TButton",
            command=lambda: self.start_execution(dry_run=True),
        )
        self.btn_dryrun.pack(side="left", padx=(0, 12))

        self.btn_run = ttk.Button(
            btn_box,
            text="[RUN] Start Pipeline",
            style="Primary.TButton",
            command=lambda: self.start_execution(dry_run=False),
        )
        self.btn_run.pack(side="left", padx=(0, 12))

        self.btn_stop = ttk.Button(
            btn_box,
            text="[STOP] Abort Execution",
            style="Danger.TButton",
            command=self.stop_execution,
            state="disabled",
        )
        self.btn_stop.pack(side="left")

        # Scrollable main configuration area
        scroll_outer = tk.Frame(parent, bg=COLOR_BG)
        scroll_outer.pack(side="top", fill="both", expand=True, padx=6, pady=(6, 2))

        canvas = tk.Canvas(scroll_outer, bg=COLOR_BG, highlightthickness=0)
        self._config_canvas = canvas
        v_scroll = ttk.Scrollbar(scroll_outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=v_scroll.set)

        v_scroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        container = ttk.Frame(canvas)
        canvas_window = canvas.create_window((0, 0), window=container, anchor="nw")

        def _on_frame_configure(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event):
            canvas.itemconfig(canvas_window, width=event.width)

        container.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)

        # Mousewheel scrolling support (routes to active scrollable tab)
        def _on_mousewheel(event):
            active_tab = self.notebook.index(self.notebook.select())
            target_cv = self._config_canvas if active_tab == 0 else (self._profile_canvas if active_tab == 1 else None)
            if target_cv is None:
                return
            if event.num == 4:
                target_cv.yview_scroll(-1, "units")
            elif event.num == 5:
                target_cv.yview_scroll(1, "units")
            elif event.delta:
                target_cv.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        canvas.bind_all("<Button-4>", _on_mousewheel)
        canvas.bind_all("<Button-5>", _on_mousewheel)

        # 2-Column Responsive Layout inside scrollable container
        left_col = ttk.Frame(container)
        left_col.pack(side="left", fill="both", expand=True, padx=(0, 6))

        right_col = ttk.Frame(container)
        right_col.pack(side="right", fill="both", expand=True, padx=(6, 0))

        # ---------------------------------------------------------------------
        # CARD 1: File & Directory Paths + Direct SRA Input
        # ---------------------------------------------------------------------
        card_paths = tk.Frame(left_col, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_paths.pack(fill="x", pady=(0, 8))

        lbl_paths_header = tk.Label(card_paths, text="Input, Output & SRA Libraries", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_paths_header, "head", "bold")
        lbl_paths_header.pack(anchor="w", padx=14, pady=(12, 8))

        # Raw FASTQ Input directory (--input)
        self.var_input = tk.StringVar(value=str(PROJECT_DIR / "data"))
        self.create_path_row(card_paths, "Input Directory (--input):", self.var_input, is_folder=True)

        # Analysis Results Output directory (--output)
        self.var_output = tk.StringVar(value=str(PROJECT_DIR / "results"))
        self.create_path_row(card_paths, "Output Directory (--output):", self.var_output, is_folder=True)

        # Temporary scratch directory (--temp-dir)
        self.var_temp = tk.StringVar(value="/tmp")
        self.create_path_row(card_paths, "Temporary Scratch Directory (--temp-dir):", self.var_temp, is_folder=True)

        # Direct SRA Accessions Text Box (New Feature: type or paste SRA IDs directly!)
        sra_direct_frame = tk.Frame(card_paths, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        sra_direct_frame.pack(fill="x", padx=14, pady=(8, 6))

        sra_top_row = tk.Frame(sra_direct_frame, bg="#172033")
        sra_top_row.pack(fill="x", padx=10, pady=(8, 2))

        lbl_sra_direct = tk.Label(
            sra_top_row,
            text="Direct SRA Accessions (type or paste IDs, e.g. SRR10677983, SRR25008973):",
            fg=COLOR_PRIMARY,
            bg="#172033",
        )
        self.register_scalable(lbl_sra_direct, "body", "bold")
        lbl_sra_direct.pack(side="left")

        self.sra_count_var = tk.StringVar(value="0 IDs typed")
        lbl_sra_count = tk.Label(sra_top_row, textvariable=self.sra_count_var, fg=COLOR_ACCENT, bg="#172033")
        self.register_scalable(lbl_sra_count, "small", "bold")
        lbl_sra_count.pack(side="right")

        sra_txt_box = tk.Frame(sra_direct_frame, bg="#172033")
        sra_txt_box.pack(fill="x", padx=10, pady=(2, 6))

        self.txt_sra_direct = tk.Text(
            sra_txt_box,
            height=3,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
            wrap="word",
        )
        self.register_scalable(self.txt_sra_direct, "mono")
        self.txt_sra_direct.pack(side="left", fill="x", expand=True, padx=(0, 8))
        self.txt_sra_direct.bind("<KeyRelease>", self._update_sra_direct_count)

        sra_btn_col = tk.Frame(sra_txt_box, bg="#172033")
        sra_btn_col.pack(side="right", fill="y")

        btn_load_sra_txt = ttk.Button(
            sra_btn_col,
            text="Load .txt...",
            style="Secondary.TButton",
            command=self.import_sra_file_into_textbox,
        )
        btn_load_sra_txt.pack(fill="x", pady=(0, 4))

        btn_clear_sra = ttk.Button(
            sra_btn_col,
            text="Clear IDs",
            style="Secondary.TButton",
            command=self.clear_sra_direct_textbox,
        )
        btn_clear_sra.pack(fill="x")

        # Optional SRA Accessions file list (--sra)
        self.var_sra = tk.StringVar(value="")
        self.create_path_row(card_paths, "Or SRA Accessions File (--sra, optional .txt):", self.var_sra, is_folder=False)

        # ---------------------------------------------------------------------
        # CARD 2: External Reference Databases
        # ---------------------------------------------------------------------
        card_dbs = tk.Frame(left_col, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_dbs.pack(fill="x", pady=6)

        lbl_dbs_header = tk.Label(card_dbs, text="External Databases (Optional - Auto-Downloads Base DBs if Empty)", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_dbs_header, "head", "bold")
        lbl_dbs_header.pack(anchor="w", padx=14, pady=(12, 4))

        lbl_dbs_sub = tk.Label(
            card_dbs,
            text="If left blank, DeepVir automatically downloads the base tutorial DBs on first run (E. coli K-12 DIAMOND & RefSeq Viral Kraken2):",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_dbs_sub, "small")
        lbl_dbs_sub.pack(anchor="w", padx=14, pady=(0, 6))

        # DIAMOND binary database (.dmnd)
        self.var_diamond_db = tk.StringVar(value="")
        self.create_path_row(
            card_dbs,
            "DIAMOND Database (.dmnd):",
            self.var_diamond_db,
            is_folder=False,
            filetypes=[("DIAMOND DB", "*.dmnd"), ("All Files", "*.*")],
        )

        # Kraken2 taxonomic database directory
        self.var_kraken2_db = tk.StringVar(value="")
        self.create_path_row(card_dbs, "Kraken2 Database (Folder):", self.var_kraken2_db, is_folder=True)

        # ---------------------------------------------------------------------
        # CARD 3: Execution Resources & Profile
        # ---------------------------------------------------------------------
        card_exec = tk.Frame(left_col, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_exec.pack(fill="x", pady=6)

        lbl_exec_header = tk.Label(card_exec, text="Resources & Snakemake Profile", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_exec_header, "head", "bold")
        lbl_exec_header.pack(anchor="w", padx=14, pady=(12, 8))

        row_params = tk.Frame(card_exec, bg=COLOR_CARD)
        row_params.pack(fill="x", padx=14, pady=(0, 10))

        # Snakemake profile selector
        lbl_prof = tk.Label(row_params, text="Profile:", bg=COLOR_CARD, fg=COLOR_TEXT)
        self.register_scalable(lbl_prof, "body")
        lbl_prof.pack(side="left", padx=(0, 8))

        self.var_profile = tk.StringVar(value="local")
        self.cb_profile = ttk.Combobox(
            row_params,
            textvariable=self.var_profile,
            values=self._get_available_profiles(),
            width=16,
            state="readonly",
        )
        self.register_scalable(self.cb_profile, "body")
        self.cb_profile.pack(side="left", padx=(0, 20))
        self.cb_profile.bind("<<ComboboxSelected>>", self._on_profile_combobox_changed)

        # Maximum concurrent jobs / threads
        lbl_jobs = tk.Label(row_params, text="Max Jobs / Cores:", bg=COLOR_CARD, fg=COLOR_TEXT)
        self.register_scalable(lbl_jobs, "body")
        lbl_jobs.pack(side="left", padx=(0, 8))

        self.var_jobs = tk.IntVar(value=8)
        self.sp_jobs = tk.Spinbox(
            row_params,
            from_=1,
            to=128,
            textvariable=self.var_jobs,
            width=6,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(self.sp_jobs, "body")
        self.sp_jobs.pack(side="left")

        # Shortcut button to open Tab 2 (Per-Rule Resources Editor)
        row_prof_btn = tk.Frame(card_exec, bg=COLOR_CARD)
        row_prof_btn.pack(fill="x", padx=14, pady=(0, 12))

        btn_edit_rules = ttk.Button(
            row_prof_btn,
            text="[Configure Per-Rule Threads & Memory (Tab 2)...]",
            style="Secondary.TButton",
            command=lambda: self.notebook.select(self.tab_profile),
        )
        btn_edit_rules.pack(fill="x")

        # ---------------------------------------------------------------------
        # CARD 4: Read QC & Quality Filtering (fastp)
        # ---------------------------------------------------------------------
        card_fastp = tk.Frame(left_col, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_fastp.pack(fill="x", pady=6)

        lbl_fastp_header = tk.Label(card_fastp, text="Read QC & Quality Filtering (fastp)", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_fastp_header, "head", "bold")
        lbl_fastp_header.pack(anchor="w", padx=14, pady=(12, 4))

        lbl_fastp_sub = tk.Label(
            card_fastp,
            text="Specify length trimming thresholds and Phred quality score for fastp:",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_fastp_sub, "small")
        lbl_fastp_sub.pack(anchor="w", padx=14, pady=(0, 10))

        row_qc_params = tk.Frame(card_fastp, bg=COLOR_CARD)
        row_qc_params.pack(fill="x", padx=14, pady=(0, 10))

        # Minimum Read Length (--fastp-length)
        lbl_len = tk.Label(row_qc_params, text="Min Length (bp):", bg=COLOR_CARD, fg=COLOR_TEXT)
        self.register_scalable(lbl_len, "body")
        lbl_len.pack(side="left", padx=(0, 8))

        self.var_fastp_len = tk.IntVar(value=30)
        self.sp_fastp_len = tk.Spinbox(
            row_qc_params,
            from_=15,
            to=50000,
            increment=5,
            textvariable=self.var_fastp_len,
            width=8,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(self.sp_fastp_len, "body")
        self.sp_fastp_len.pack(side="left", padx=(0, 24))

        # Minimum Quality Phred (--fastp-quality)
        lbl_qual = tk.Label(row_qc_params, text="Min Phred Quality (Q):", bg=COLOR_CARD, fg=COLOR_TEXT)
        self.register_scalable(lbl_qual, "body")
        lbl_qual.pack(side="left", padx=(0, 8))

        self.var_fastp_qual = tk.IntVar(value=30)
        self.sp_fastp_qual = tk.Spinbox(
            row_qc_params,
            from_=1,
            to=45,
            increment=1,
            textvariable=self.var_fastp_qual,
            width=6,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(self.sp_fastp_qual, "body")
        self.sp_fastp_qual.pack(side="left")

        # Quick Presets Buttons row
        row_presets = tk.Frame(card_fastp, bg=COLOR_CARD)
        row_presets.pack(fill="x", padx=14, pady=(0, 14))

        lbl_presets = tk.Label(row_presets, text="Presets:", bg=COLOR_CARD, fg=COLOR_MUTED)
        self.register_scalable(lbl_presets, "small", "bold")
        lbl_presets.pack(side="left", padx=(0, 8))

        def set_fastp_preset(length, qual):
            self.var_fastp_len.set(length)
            self.var_fastp_qual.set(qual)

        btn_illumina = tk.Button(
            row_presets,
            text="Illumina (L: 30, Q: 30)",
            bg="#334155",
            fg=COLOR_TEXT,
            relief="flat",
            padx=10,
            pady=4,
            cursor="hand2",
            activebackground="#475569",
            activeforeground="#ffffff",
            command=lambda: set_fastp_preset(30, 30),
        )
        self.register_scalable(btn_illumina, "small", "bold")
        btn_illumina.pack(side="left", padx=(0, 8))

        btn_nanopore = tk.Button(
            row_presets,
            text="Nanopore / Long-read (L: 500, Q: 10)",
            bg="#334155",
            fg=COLOR_TEXT,
            relief="flat",
            padx=10,
            pady=4,
            cursor="hand2",
            activebackground="#475569",
            activeforeground="#ffffff",
            command=lambda: set_fastp_preset(500, 10),
        )
        self.register_scalable(btn_nanopore, "small", "bold")
        btn_nanopore.pack(side="left")

        # ---------------------------------------------------------------------
        # RIGHT COLUMN: Analysis Modules Selection Card
        # ---------------------------------------------------------------------
        card_modules = tk.Frame(right_col, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_modules.pack(fill="both", expand=True, pady=(0, 8))

        lbl_mod_header = tk.Label(card_modules, text="Analysis Modules Selection", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_mod_header, "head", "bold")
        lbl_mod_header.pack(anchor="w", padx=16, pady=(12, 4))

        lbl_mod_sub = tk.Label(
            card_modules,
            text="Select pipeline analysis modules to enable:",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_mod_sub, "small")
        lbl_mod_sub.pack(anchor="w", padx=16, pady=(0, 10))

        # Checkbox state variables
        self.mod_host_filter = tk.BooleanVar(value=False)
        self.var_host_accession = tk.StringVar(value="GCA_022682495.2, GCF_000819615.1")
        self.mod_reads_kraken = tk.BooleanVar(value=False)
        self.mod_reads_diamond = tk.BooleanVar(value=False)
        self.mod_assembly = tk.BooleanVar(value=True)
        self.mod_kraken2 = tk.BooleanVar(value=True)
        self.mod_diamond = tk.BooleanVar(value=True)
        self.mod_darkmatter = tk.BooleanVar(value=True)
        self.var_rdrp_minscore = tk.IntVar(value=75)
        self.mod_rvdb = tk.BooleanVar(value=False)
        self.mod_download_only = tk.BooleanVar(value=False)
        self.mod_remove_download = tk.BooleanVar(value=False)

        # De novo Assembler engine selection (config.yaml -> tool: denovo)
        self.asm_spades = tk.BooleanVar(value=True)
        self.asm_megahit = tk.BooleanVar(value=False)
        self.asm_flye = tk.BooleanVar(value=False)
        self.asm_raven = tk.BooleanVar(value=False)
        self.asm_medaka_flye = tk.BooleanVar(value=False)
        self.asm_medaka_raven = tk.BooleanVar(value=False)

        # Host Genome Filtering / Depletion Group (BWA)
        self.build_host_filter_group(card_modules)

        # Categorized module groups
        self.create_check_group(
            card_modules,
            "> Read-Level Taxonomy (Fast Screening)",
            [
                ("Taxonomic classification of reads via Kraken2 (--reads-kraken)", self.mod_reads_kraken),
                ("Reads alignment via DIAMOND + MEGAN LCA (--reads-diamond)", self.mod_reads_diamond),
            ],
        )

        # De Novo Assembly & Contigs Group (with nested assembler tools selector)
        self.build_assembly_group(card_modules)

        # Viral Discovery / Dark Matter Group (with configurable RdRp minscore)
        self.build_darkmatter_group(card_modules)

        self.create_check_group(
            card_modules,
            "> SRA Data Management",
            [
                ("Download SRA accessions only without running analysis (--download-only)", self.mod_download_only),
                ("Remove downloaded raw SRA archives after extraction (--remove-download)", self.mod_remove_download),
            ],
        )

    # =========================================================================
    # TAB 2: PROFILE & PER-RULE RESOURCES BUILDER
    # =========================================================================
    def build_tab_profile(self, parent):
        """
        Builds the Profile & Rule Resources tab:
        - Allows selecting the active profile (`local` vs `profile_slurm`).
        - Edits default resources (`threads`, `mem_mb`, `restart-times`, `runtime`, `latency-wait`).
        - Provides 1-click Hardware Presets (Laptop, Desktop, Server, Slurm HPC).
        - Displays an interactive grid of per-rule `threads` and `mem_mb` (`set-resources`).
        - Saves directly to `pipeline/profiles/<profile>/config.yaml`.
        """
        # Bottom Action Bar for Saving / Reloading Profile
        bottom_bar = tk.Frame(parent, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        bottom_bar.pack(side="bottom", fill="x", padx=6, pady=(4, 6))

        bot_inner = tk.Frame(bottom_bar, bg=COLOR_CARD)
        bot_inner.pack(fill="x", padx=16, pady=10)

        btn_save_prof = ttk.Button(
            bot_inner,
            text="[SAVE] Save Resources to Profile (config.yaml)",
            style="Primary.TButton",
            command=lambda: self.save_editor_to_profile(show_popup=True),
        )
        btn_save_prof.pack(side="left", padx=(0, 12))

        btn_reload_prof = ttk.Button(
            bot_inner,
            text="[Reload from Disk]",
            style="Secondary.TButton",
            command=lambda: self.load_profile_into_editor(self.var_profile.get(), silent=False),
        )
        btn_reload_prof.pack(side="left", padx=(0, 12))

        self.lbl_prof_path = tk.Label(
            bot_inner,
            text="File: pipeline/profiles/local/config.yaml",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(self.lbl_prof_path, "small")
        self.lbl_prof_path.pack(side="right")

        # Scrollable main container for Profile & Rule Resources
        scroll_outer = tk.Frame(parent, bg=COLOR_BG)
        scroll_outer.pack(side="top", fill="both", expand=True, padx=6, pady=(6, 2))

        canvas = tk.Canvas(scroll_outer, bg=COLOR_BG, highlightthickness=0)
        self._profile_canvas = canvas
        v_scroll = ttk.Scrollbar(scroll_outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=v_scroll.set)

        v_scroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        container = ttk.Frame(canvas)
        canvas_window = canvas.create_window((0, 0), window=container, anchor="nw")

        def _on_frame_configure(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event):
            canvas.itemconfig(canvas_window, width=event.width)

        container.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)

        # Top Card: Profile Selection, Global Defaults & Hardware Presets
        card_top = tk.Frame(container, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_top.pack(fill="x", pady=(0, 8))

        lbl_top_title = tk.Label(
            card_top,
            text="Snakemake Profile & Default Resource Limits",
            fg=COLOR_ACCENT,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_top_title, "head", "bold")
        lbl_top_title.pack(anchor="w", padx=16, pady=(12, 4))

        lbl_top_sub = tk.Label(
            card_top,
            text="Select the target profile and configure global fallback resources and per-rule CPU/RAM allocations:",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_top_sub, "small")
        lbl_top_sub.pack(anchor="w", padx=16, pady=(0, 10))

        # Row 1: Profile selector + Default Threads + Default Mem + Restart times + Slurm Runtime
        row_globals = tk.Frame(card_top, bg=COLOR_CARD)
        row_globals.pack(fill="x", padx=16, pady=(0, 10))

        lbl_p = tk.Label(row_globals, text="Active Profile:", fg=COLOR_PRIMARY, bg=COLOR_CARD)
        self.register_scalable(lbl_p, "body", "bold")
        lbl_p.pack(side="left", padx=(0, 6))

        self.cb_profile_tab2 = ttk.Combobox(
            row_globals,
            textvariable=self.var_profile,
            values=self._get_available_profiles(),
            width=15,
            state="readonly",
        )
        self.register_scalable(self.cb_profile_tab2, "body")
        self.cb_profile_tab2.pack(side="left", padx=(0, 18))
        self.cb_profile_tab2.bind("<<ComboboxSelected>>", self._on_profile_combobox_changed)

        lbl_dt = tk.Label(row_globals, text="Default Threads:", fg=COLOR_TEXT, bg=COLOR_CARD)
        self.register_scalable(lbl_dt, "body")
        lbl_dt.pack(side="left", padx=(0, 6))

        self.var_def_threads = tk.IntVar(value=1)
        sp_dt = tk.Spinbox(
            row_globals,
            from_=1,
            to=64,
            textvariable=self.var_def_threads,
            width=5,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(sp_dt, "body")
        sp_dt.pack(side="left", padx=(0, 18))

        lbl_dm = tk.Label(row_globals, text="Default Memory (MB):", fg=COLOR_TEXT, bg=COLOR_CARD)
        self.register_scalable(lbl_dm, "body")
        lbl_dm.pack(side="left", padx=(0, 6))

        self.var_def_mem = tk.IntVar(value=2000)
        sp_dm = tk.Spinbox(
            row_globals,
            from_=500,
            to=128000,
            increment=500,
            textvariable=self.var_def_mem,
            width=7,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(sp_dm, "body")
        sp_dm.pack(side="left", padx=(0, 18))

        lbl_rt = tk.Label(row_globals, text="Retries:", fg=COLOR_TEXT, bg=COLOR_CARD)
        self.register_scalable(lbl_rt, "body")
        lbl_rt.pack(side="left", padx=(0, 6))

        self.var_restarts = tk.IntVar(value=2)
        sp_rt = tk.Spinbox(
            row_globals,
            from_=0,
            to=10,
            textvariable=self.var_restarts,
            width=4,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(sp_rt, "body")
        sp_rt.pack(side="left", padx=(0, 18))

        lbl_wall = tk.Label(row_globals, text="Slurm Walltime:", fg=COLOR_TEXT, bg=COLOR_CARD)
        self.register_scalable(lbl_wall, "body")
        lbl_wall.pack(side="left", padx=(0, 6))

        self.var_slurm_runtime = tk.StringVar(value="1-00:00:00")
        ent_wall = tk.Entry(
            row_globals,
            textvariable=self.var_slurm_runtime,
            width=11,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )
        self.register_scalable(ent_wall, "body")
        ent_wall.pack(side="left")

        # Row 2: 1-Click Hardware Presets
        row_hw_presets = tk.Frame(card_top, bg=COLOR_CARD)
        row_hw_presets.pack(fill="x", padx=16, pady=(2, 12))

        lbl_hw = tk.Label(row_hw_presets, text="Hardware Presets (Apply to All Rules):", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(lbl_hw, "small", "bold")
        lbl_hw.pack(side="left", padx=(0, 10))

        for label_txt, threads_val, mem_val, def_mem in [
            ("Laptop (4 threads / 8 GB)", 4, 8000, 2000),
            ("Desktop Default (6 threads / 12 GB)", 6, 12000, 2000),
            ("Workstation (12 threads / 32 GB)", 12, 32000, 4000),
            ("HPC Slurm (16 threads / 60 GB)", 16, 60000, 4000),
        ]:
            btn_p = tk.Button(
                row_hw_presets,
                text=label_txt,
                bg="#334155",
                fg=COLOR_TEXT,
                relief="flat",
                padx=10,
                pady=4,
                cursor="hand2",
                activebackground="#475569",
                activeforeground="#ffffff",
                command=lambda t=threads_val, m=mem_val, dm=def_mem: self.apply_hardware_preset(t, m, dm),
            )
            self.register_scalable(btn_p, "small", "bold")
            btn_p.pack(side="left", padx=4)

        # Bottom Card: Per-Rule Resources Table (`set-resources`)
        card_rules = tk.Frame(container, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card_rules.pack(fill="both", expand=True, pady=(0, 8))

        lbl_rules_title = tk.Label(
            card_rules,
            text="Per-Rule Resource Allocation (set-resources)",
            fg=COLOR_ACCENT,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_rules_title, "head", "bold")
        lbl_rules_title.pack(anchor="w", padx=16, pady=(12, 8))

        # Table Header
        hdr_row = tk.Frame(card_rules, bg="#0f172a", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        hdr_row.pack(fill="x", padx=16, pady=(0, 4))

        lbl_h1 = tk.Label(hdr_row, text="Snakemake Rule", width=26, anchor="w", fg=COLOR_PRIMARY, bg="#0f172a")
        self.register_scalable(lbl_h1, "body", "bold")
        lbl_h1.pack(side="left", padx=10, pady=6)

        lbl_h2 = tk.Label(hdr_row, text="Module / Description", width=36, anchor="w", fg=COLOR_PRIMARY, bg="#0f172a")
        self.register_scalable(lbl_h2, "body", "bold")
        lbl_h2.pack(side="left", padx=10, pady=6)

        lbl_h3 = tk.Label(hdr_row, text="Threads / CPUs", width=15, anchor="center", fg=COLOR_PRIMARY, bg="#0f172a")
        self.register_scalable(lbl_h3, "body", "bold")
        lbl_h3.pack(side="left", padx=10, pady=6)

        lbl_h4 = tk.Label(hdr_row, text="Memory (MB / GB)", width=24, anchor="w", fg=COLOR_PRIMARY, bg="#0f172a")
        self.register_scalable(lbl_h4, "body", "bold")
        lbl_h4.pack(side="left", padx=10, pady=6)

        # Rows for each rule in RULE_DEFINITIONS
        for idx, (rule_name, rule_desc, category) in enumerate(RULE_DEFINITIONS):
            row_bg = "#172033" if idx % 2 == 0 else COLOR_CARD
            r_frame = tk.Frame(card_rules, bg=row_bg)
            r_frame.pack(fill="x", padx=16, pady=1)

            lbl_rname = tk.Label(r_frame, text=rule_name, width=26, anchor="w", fg=COLOR_TEXT, bg=row_bg)
            self.register_scalable(lbl_rname, "mono", "bold")
            lbl_rname.pack(side="left", padx=10, pady=5)

            lbl_rdesc = tk.Label(r_frame, text=f"[{category}] {rule_desc}", width=36, anchor="w", fg=COLOR_MUTED, bg=row_bg)
            self.register_scalable(lbl_rdesc, "small")
            lbl_rdesc.pack(side="left", padx=10, pady=5)

            # Threads Spinbox
            t_box = tk.Frame(r_frame, bg=row_bg, width=140)
            t_box.pack(side="left", padx=10, pady=3)
            var_t = tk.IntVar(value=6)
            sp_t = tk.Spinbox(
                t_box,
                from_=1,
                to=128,
                textvariable=var_t,
                width=6,
                bg="#0f172a",
                fg=COLOR_TEXT,
                insertbackground="white",
            )
            self.register_scalable(sp_t, "body")
            sp_t.pack(side="left", padx=12)

            # Memory (mem_mb) Spinbox + live GB indicator
            m_box = tk.Frame(r_frame, bg=row_bg)
            m_box.pack(side="left", padx=10, pady=3)
            var_m = tk.IntVar(value=12000)
            lbl_gb = tk.Label(m_box, text="(12.0 GB)", fg=COLOR_ACCENT, bg=row_bg, width=10, anchor="w")
            self.register_scalable(lbl_gb, "small", "bold")

            def _make_mem_trace(v_m=var_m, l_gb=lbl_gb):
                def _on_mem_change(*_args):
                    try:
                        mb = int(v_m.get())
                        l_gb.config(text=f"({mb / 1000.0:.1f} GB)")
                    except (ValueError, tk.TclError):
                        pass
                return _on_mem_change

            var_m.trace_add("write", _make_mem_trace())

            sp_m = tk.Spinbox(
                m_box,
                from_=1000,
                to=512000,
                increment=1000,
                textvariable=var_m,
                width=8,
                bg="#0f172a",
                fg=COLOR_TEXT,
                insertbackground="white",
            )
            self.register_scalable(sp_m, "body")
            sp_m.pack(side="left", padx=(0, 6))
            lbl_gb.pack(side="left")

            self.rule_resource_vars[rule_name] = (var_t, var_m, lbl_gb)

        tk.Frame(card_rules, bg=COLOR_CARD, height=10).pack()

    def _get_available_profiles(self) -> list[str]:
        """Returns available profile folder names inside `pipeline/profiles/`."""
        profiles = ["local", "profile_slurm"]
        if PROFILES_DIR.is_dir():
            for d in sorted(PROFILES_DIR.iterdir()):
                if d.is_dir() and d.name not in profiles:
                    profiles.append(d.name)
        return profiles

    def _on_profile_combobox_changed(self, _event=None):
        """Synchronizes mode radio buttons and loads the selected profile into Tab 2."""
        prof = self.var_profile.get().strip()
        if "slurm" in prof:
            self.mode_var.set("slurm")
        else:
            self.mode_var.set("local")
        self.load_profile_into_editor(prof, silent=True)

    def load_profile_into_editor(self, profile_name: str, silent: bool = False):
        """Loads `pipeline/profiles/<profile_name>/config.yaml` into Tab 2's spinboxes."""
        prof_path = PROFILES_DIR / profile_name / "config.yaml"
        data = parse_profile_yaml(prof_path)

        self.var_jobs.set(data.get("jobs", 8))
        self.var_def_threads.set(data.get("default_threads", 1))
        self.var_def_mem.set(data.get("default_mem_mb", 2000))
        self.var_restarts.set(data.get("restart_times", 2))
        self.var_slurm_runtime.set(data.get("runtime", "1-00:00:00"))

        if hasattr(self, "lbl_prof_path"):
            self.lbl_prof_path.config(text=f"File: pipeline/profiles/{profile_name}/config.yaml")

        rules_data = data.get("rules", {})
        for rule_name, (var_t, var_m, _lbl_gb) in self.rule_resource_vars.items():
            rcfg = rules_data.get(rule_name, {"threads": 6, "mem_mb": 12000})
            var_t.set(rcfg.get("threads", 6))
            var_m.set(rcfg.get("mem_mb", 12000))

        if not silent:
            self.append_terminal(f"[PROFILE LOADED] Loaded rule resources from pipeline/profiles/{profile_name}/config.yaml\n")

    def save_editor_to_profile(self, show_popup: bool = False):
        """Writes the current Tab 2 resource values back to `pipeline/profiles/<profile>/config.yaml`."""
        profile_name = self.var_profile.get().strip() or "local"
        prof_path = PROFILES_DIR / profile_name / "config.yaml"

        rules_dict = {}
        for rule_name, (var_t, var_m, _) in self.rule_resource_vars.items():
            try:
                t_val = max(1, int(var_t.get()))
            except (ValueError, tk.TclError):
                t_val = 6
            try:
                m_val = max(500, int(var_m.get()))
            except (ValueError, tk.TclError):
                m_val = 12000
            rules_dict[rule_name] = {"threads": t_val, "mem_mb": m_val}

        try:
            jobs_val = max(1, int(self.var_jobs.get()))
        except (ValueError, tk.TclError):
            jobs_val = 8

        try:
            def_t = max(1, int(self.var_def_threads.get()))
        except (ValueError, tk.TclError):
            def_t = 1

        try:
            def_m = max(500, int(self.var_def_mem.get()))
        except (ValueError, tk.TclError):
            def_m = 2000

        try:
            restarts = max(0, int(self.var_restarts.get()))
        except (ValueError, tk.TclError):
            restarts = 2

        data = {
            "jobs": jobs_val,
            "restart_times": restarts,
            "latency_wait": 600,
            "runtime": self.var_slurm_runtime.get().strip() or "1-00:00:00",
            "is_slurm": "slurm" in profile_name.lower(),
            "default_threads": def_t,
            "default_mem_mb": def_m,
            "rules": rules_dict,
        }

        write_profile_yaml(prof_path, data)
        self.append_terminal(f"[PROFILE SAVED] Updated per-rule resources in {prof_path}\n")
        if show_popup:
            messagebox.showinfo(
                "Profile Saved",
                f"Per-rule resources and profile settings were saved to:\n{prof_path}",
            )

    def apply_hardware_preset(self, threads_val: int, mem_val: int, def_mem: int):
        """Applies a 1-click hardware preset to all rule spinboxes in Tab 2."""
        self.var_def_mem.set(def_mem)
        for rule_name, (var_t, var_m, _) in self.rule_resource_vars.items():
            if rule_name == "cd_hit":
                var_t.set(min(4, threads_val))
                var_m.set(min(16000, mem_val))
            else:
                var_t.set(threads_val)
                var_m.set(mem_val)

    # =========================================================================
    # DIRECT SRA ACCESSIONS HELPERS
    # =========================================================================
    def get_direct_sra_tokens(self) -> list[str]:
        """Extracts valid SRA accession IDs typed directly in the GUI text box."""
        if not hasattr(self, "txt_sra_direct"):
            return []
        raw = self.txt_sra_direct.get("1.0", "end-1c").strip()
        if not raw:
            return []
        tokens = []
        for tok in re.split(r"[,;\s]+", raw):
            clean = tok.strip()
            if clean and not clean.startswith("#") and clean not in tokens:
                tokens.append(clean)
        return tokens

    def _update_sra_direct_count(self, _event=None):
        tokens = self.get_direct_sra_tokens()
        self.sra_count_var.set(f"{len(tokens)} ID{'s' if len(tokens) != 1 else ''} typed")

    def clear_sra_direct_textbox(self):
        self.txt_sra_direct.delete("1.0", "end")
        self._update_sra_direct_count()

    def import_sra_file_into_textbox(self):
        """Loads an external `.txt` file of SRA accessions directly into the GUI text box."""
        chosen = filedialog.askopenfilename(
            initialdir=str(PROJECT_DIR / "data"),
            title="Select SRA Accessions File (.txt)",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
        )
        if chosen and Path(chosen).is_file():
            content = Path(chosen).read_text(encoding="utf-8", errors="replace").strip()
            existing = self.txt_sra_direct.get("1.0", "end-1c").strip()
            combined = f"{existing}\n{content}".strip() if existing else content
            self.txt_sra_direct.delete("1.0", "end")
            self.txt_sra_direct.insert("1.0", combined)
            self._update_sra_direct_count()

    def create_path_row(self, parent, label_text, str_var, is_folder=True, filetypes=None):
        """
        Creates an interactive path row with enlarged label, entry field, and Browse button.
        """
        frame = tk.Frame(parent, bg=COLOR_CARD)
        frame.pack(fill="x", padx=14, pady=5)

        lbl = tk.Label(frame, text=label_text, fg=COLOR_TEXT, bg=COLOR_CARD)
        self.register_scalable(lbl, "body", "bold")
        lbl.pack(anchor="w")

        input_box = tk.Frame(frame, bg=COLOR_CARD)
        input_box.pack(fill="x", pady=(3, 4))

        entry = tk.Entry(
            input_box,
            textvariable=str_var,
            bg="#0f172a",
            fg=COLOR_TEXT,
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )
        self.register_scalable(entry, "body")
        entry.pack(side="left", fill="x", expand=True, ipady=6, padx=(0, 8))

        def browse():
            if is_folder:
                chosen = filedialog.askdirectory(initialdir=str_var.get() or str(PROJECT_DIR))
            else:
                chosen = filedialog.askopenfilename(
                    initialdir=str_var.get() or str(PROJECT_DIR),
                    filetypes=filetypes or [("All Files", "*.*")],
                )
            if chosen:
                str_var.set(chosen)

        btn_browse = ttk.Button(input_box, text="Browse...", style="Secondary.TButton", command=browse)
        btn_browse.pack(side="right")

    def create_check_group(self, parent, title, items):
        """
        Renders a grouped box of checkboxes with clean ASCII titles and enlarged fonts.
        """
        box = tk.Frame(parent, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        box.pack(fill="x", padx=14, pady=5)

        lbl = tk.Label(box, text=title, fg=COLOR_ACCENT, bg="#172033")
        self.register_scalable(lbl, "subhead", "bold")
        lbl.pack(anchor="w", padx=10, pady=(8, 4))

        for text, var in items:
            chk = ttk.Checkbutton(box, text=text, variable=var)
            chk.pack(anchor="w", padx=14, pady=3)

        tk.Frame(box, bg="#172033", height=6).pack()

    def build_host_filter_group(self, parent):
        """
        Constructs the Host Genome Filtering / Depletion (BWA) boolean module section.
        Supports one or multiple target NCBI host genome accessions (comma-separated).
        Downloads each specified NCBI host assembly, indexes them together with BWA,
        quantifies total & per-host read percentages, and removes host reads before downstream steps.
        """
        box = tk.Frame(parent, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        box.pack(fill="x", padx=14, pady=5)

        lbl = tk.Label(box, text="> Host Genome Filtering / Depletion (BWA - Multi-Host Ready)", fg=COLOR_ACCENT, bg="#172033")
        self.register_scalable(lbl, "subhead", "bold")
        lbl.pack(anchor="w", padx=10, pady=(8, 4))

        chk_hf = ttk.Checkbutton(
            box,
            text="Filter out host reads via BWA & quantify host % (--host-filter)",
            variable=self.mod_host_filter,
            command=self.on_host_filter_toggle,
        )
        chk_hf.pack(anchor="w", padx=14, pady=3)

        self.frame_host_acc = tk.Frame(box, bg="#0f172a", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        self.frame_host_acc.pack(fill="x", padx=16, pady=(4, 8))

        acc_top_row = tk.Frame(self.frame_host_acc, bg="#0f172a")
        acc_top_row.pack(fill="x", padx=10, pady=(6, 2))

        lbl_acc = tk.Label(
            acc_top_row,
            text="Target Host NCBI Accession(s) (1 or more, comma-separated):",
            fg=COLOR_PRIMARY,
            bg="#0f172a",
        )
        self.register_scalable(lbl_acc, "body", "bold")
        lbl_acc.pack(side="left")

        self.host_count_var = tk.StringVar(value="1 target genome")
        lbl_host_cnt = tk.Label(acc_top_row, textvariable=self.host_count_var, fg=COLOR_ACCENT, bg="#0f172a")
        self.register_scalable(lbl_host_cnt, "small", "bold")
        lbl_host_cnt.pack(side="right")

        acc_row = tk.Frame(self.frame_host_acc, bg="#0f172a")
        acc_row.pack(fill="x", padx=10, pady=(2, 4))

        self.ent_host_acc = tk.Entry(
            acc_row,
            textvariable=self.var_host_accession,
            bg="#1e293b",
            fg=COLOR_TEXT,
            disabledbackground="#0f172a",
            disabledforeground="#64748b",
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )
        self.register_scalable(self.ent_host_acc, "mono", "bold")
        self.ent_host_acc.pack(side="left", fill="x", expand=True, ipady=4)
        self.var_host_accession.trace_add("write", lambda *_: self._update_host_acc_count())

        preset_row1 = tk.Frame(self.frame_host_acc, bg="#0f172a")
        preset_row1.pack(fill="x", padx=10, pady=(2, 2))

        lbl_hp = tk.Label(preset_row1, text="Recommended Hosts:", fg=COLOR_MUTED, bg="#0f172a")
        self.register_scalable(lbl_hp, "small", "bold")
        lbl_hp.pack(side="left", padx=(0, 6))

        preset_row2 = tk.Frame(self.frame_host_acc, bg="#0f172a")
        preset_row2.pack(fill="x", padx=10, pady=(2, 6))

        self.host_preset_buttons = []
        host_presets = [
            (preset_row1, "+ PhiX174 (GCF_000819615.1)", "GCF_000819615.1"),
            (preset_row1, "+ D. rotundus (GCA_022682495.2)", "GCA_022682495.2"),
            (preset_row2, "+ H. sapiens (GCF_000001405.40)", "GCF_000001405.40"),
            (preset_row2, "+ A. aegypti (GCF_002204515.2)", "GCF_002204515.2"),
        ]
        for row_widget, label_txt, acc_code in host_presets:
            btn_hp = tk.Button(
                row_widget,
                text=label_txt,
                bg="#334155",
                fg=COLOR_TEXT,
                relief="flat",
                padx=8,
                pady=2,
                cursor="hand2",
                activebackground="#475569",
                activeforeground="#ffffff",
                command=lambda a=acc_code: self._toggle_host_accession_preset(a),
            )
            self.register_scalable(btn_hp, "small", "bold")
            btn_hp.pack(side="left", padx=3)
            self.host_preset_buttons.append(btn_hp)

        btn_clear_hp = tk.Button(
            preset_row2,
            text="Clear",
            bg="#334155",
            fg=COLOR_TEXT,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            activebackground="#475569",
            activeforeground="#ffffff",
            command=lambda: self.var_host_accession.set(""),
        )
        self.register_scalable(btn_clear_hp, "small", "bold")
        btn_clear_hp.pack(side="left", padx=3)
        self.host_preset_buttons.append(btn_clear_hp)

        tk.Frame(box, bg="#172033", height=4).pack()
        self._update_host_acc_count()
        self.on_host_filter_toggle()

    def get_host_accession_tokens(self) -> list[str]:
        """Returns the list of unique target host NCBI accessions entered in the GUI."""
        raw = self.var_host_accession.get().strip()
        tokens = []
        for tok in re.split(r"[,;\s]+", raw):
            clean = tok.strip()
            if clean and clean not in tokens:
                tokens.append(clean)
        return tokens

    def _update_host_acc_count(self):
        if not hasattr(self, "host_count_var"):
            return
        tokens = self.get_host_accession_tokens()
        n = len(tokens)
        self.host_count_var.set(f"{n} target genome{'s' if n != 1 else ''}")

    def _toggle_host_accession_preset(self, acc_code: str):
        """Adds `acc_code` to the comma-separated host accession field (or removes it if already present)."""
        tokens = self.get_host_accession_tokens()
        if acc_code in tokens:
            tokens.remove(acc_code)
        else:
            tokens.append(acc_code)
        self.var_host_accession.set(", ".join(tokens))

    def on_host_filter_toggle(self):
        state = "normal" if self.mod_host_filter.get() else "disabled"
        self.ent_host_acc.config(state=state)
        for btn in getattr(self, "host_preset_buttons", []):
            btn.config(state=state)

    def build_assembly_group(self, parent):
        """
        Constructs the De Novo Assembly & Contigs module section.
        """
        box = tk.Frame(parent, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        box.pack(fill="x", padx=14, pady=5)

        lbl = tk.Label(box, text="> De Novo Assembly & Contigs", fg=COLOR_ACCENT, bg="#172033")
        self.register_scalable(lbl, "subhead", "bold")
        lbl.pack(anchor="w", padx=10, pady=(8, 4))

        # Main Assembly Toggle
        chk_asm = ttk.Checkbutton(
            box,
            text="De Novo Reads Assembly (--assembly)",
            variable=self.mod_assembly,
            command=self.on_assembly_toggle,
        )
        chk_asm.pack(anchor="w", padx=14, pady=3)

        # Sub-card: Assembler selection (tool.denovo in config.yaml)
        self.frame_assemblers = tk.Frame(box, bg="#0f172a", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        self.frame_assemblers.pack(fill="x", padx=16, pady=(4, 8))

        lbl_tools = tk.Label(
            self.frame_assemblers,
            text="Assembler Engine Selection (tool.denovo):",
            fg=COLOR_PRIMARY,
            bg="#0f172a",
        )
        self.register_scalable(lbl_tools, "body", "bold")
        lbl_tools.pack(anchor="w", padx=10, pady=(6, 4))

        grid_tools = tk.Frame(self.frame_assemblers, bg="#0f172a")
        grid_tools.pack(fill="x", padx=10, pady=(0, 6))

        # Column 1: Short Reads & Hybrids
        col1 = tk.Frame(grid_tools, bg="#0f172a")
        col1.pack(side="left", fill="x", expand=True)

        self.chk_spades = ttk.Checkbutton(col1, text="SPAdes (Illumina RNA/viral - default)", variable=self.asm_spades)
        self.chk_spades.pack(anchor="w", pady=2)

        self.chk_megahit = ttk.Checkbutton(col1, text="MEGAHIT (Illumina - fast / low RAM)", variable=self.asm_megahit)
        self.chk_megahit.pack(anchor="w", pady=2)

        self.chk_medaka_flye = ttk.Checkbutton(col1, text="Medaka + Flye (ONP polished)", variable=self.asm_medaka_flye)
        self.chk_medaka_flye.pack(anchor="w", pady=2)

        # Column 2: Long Reads
        col2 = tk.Frame(grid_tools, bg="#0f172a")
        col2.pack(side="left", fill="x", expand=True)

        self.chk_flye = ttk.Checkbutton(col2, text="Flye (Nanopore / PacBio long reads)", variable=self.asm_flye)
        self.chk_flye.pack(anchor="w", pady=2)

        self.chk_raven = ttk.Checkbutton(col2, text="Raven (Nanopore long reads - fast)", variable=self.asm_raven)
        self.chk_raven.pack(anchor="w", pady=2)

        self.chk_medaka_raven = ttk.Checkbutton(col2, text="Medaka + Raven (ONP polished)", variable=self.asm_medaka_raven)
        self.chk_medaka_raven.pack(anchor="w", pady=2)

        # Downstream contig annotation toggles
        chk_k2 = ttk.Checkbutton(box, text="Taxonomic classification of contigs via Kraken2 (--kraken2)", variable=self.mod_kraken2)
        chk_k2.pack(anchor="w", padx=14, pady=3)

        chk_dia = ttk.Checkbutton(box, text="Protein annotation of contigs via DIAMOND NR (--diamond)", variable=self.mod_diamond)
        chk_dia.pack(anchor="w", padx=14, pady=3)

        tk.Frame(box, bg="#172033", height=6).pack()

        # Initial synchronization of enabled state
        self.on_assembly_toggle()

    def on_assembly_toggle(self):
        state = "normal" if self.mod_assembly.get() else "disabled"
        for chk in [
            self.chk_spades,
            self.chk_megahit,
            self.chk_flye,
            self.chk_raven,
            self.chk_medaka_flye,
            self.chk_medaka_raven,
        ]:
            chk.config(state=state)

    def build_darkmatter_group(self, parent):
        """
        Constructs the Viral Discovery / Dark Matter Module section with an interactive
        RdRp MinScore threshold selector (--rdrp-minscore, 0-100, default: 75).
        """
        box = tk.Frame(parent, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        box.pack(fill="x", padx=14, pady=5)

        lbl = tk.Label(box, text="> Viral Discovery / Dark Matter Module", fg=COLOR_ACCENT, bg="#172033")
        self.register_scalable(lbl, "subhead", "bold")
        lbl.pack(anchor="w", padx=10, pady=(8, 4))

        chk_dm = ttk.Checkbutton(
            box,
            text="RdRp catalytic motif detection via PalmAnnot (--darkmatter)",
            variable=self.mod_darkmatter,
            command=self.on_darkmatter_toggle,
        )
        chk_dm.pack(anchor="w", padx=14, pady=3)

        # Sub-card for RdRp MinScore parameter
        self.frame_dm_params = tk.Frame(box, bg="#0f172a", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        self.frame_dm_params.pack(fill="x", padx=16, pady=(4, 6))

        dm_row = tk.Frame(self.frame_dm_params, bg="#0f172a")
        dm_row.pack(fill="x", padx=10, pady=6)

        lbl_ms = tk.Label(
            dm_row,
            text="RdRp MinScore (--rdrp-minscore, 0-100):",
            fg=COLOR_PRIMARY,
            bg="#0f172a",
        )
        self.register_scalable(lbl_ms, "body", "bold")
        lbl_ms.pack(side="left", padx=(0, 8))

        self.sp_rdrp_minscore = tk.Spinbox(
            dm_row,
            from_=0,
            to=100,
            increment=5,
            textvariable=self.var_rdrp_minscore,
            width=6,
            bg="#1e293b",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(self.sp_rdrp_minscore, "mono", "bold")
        self.sp_rdrp_minscore.pack(side="left", padx=(0, 12))

        self.dm_preset_buttons = []
        for p_label, p_val in [
            ("Relaxed (50)", 50),
            ("Moderate (65)", 65),
            ("Default (75)", 75),
            ("Strict (90)", 90),
        ]:
            btn_ms = tk.Button(
                dm_row,
                text=p_label,
                bg="#334155",
                fg=COLOR_TEXT,
                relief="flat",
                padx=7,
                pady=2,
                cursor="hand2",
                activebackground="#475569",
                activeforeground="#ffffff",
                command=lambda v=p_val: self._set_rdrp_minscore_preset(v),
            )
            self.register_scalable(btn_ms, "small", "bold")
            btn_ms.pack(side="left", padx=2)
            self.dm_preset_buttons.append(btn_ms)

        self.var_rdrp_minscore.trace_add("write", lambda *_: self._sync_rdrp_minscore_to_results_tab())

        chk_rvdb = ttk.Checkbutton(
            box,
            text="Cross-validation against Reference Viral DB (--rvdb)",
            variable=self.mod_rvdb,
        )
        chk_rvdb.pack(anchor="w", padx=14, pady=3)

        tk.Frame(box, bg="#172033", height=6).pack()
        self.on_darkmatter_toggle()

    def _set_rdrp_minscore_preset(self, val: int):
        self.var_rdrp_minscore.set(val)
        self._sync_rdrp_minscore_to_results_tab()

    def _sync_rdrp_minscore_to_results_tab(self):
        if hasattr(self, "var_dm_filter_minscore"):
            try:
                self.var_dm_filter_minscore.set(int(self.var_rdrp_minscore.get()))
            except (ValueError, tk.TclError):
                pass

    def on_darkmatter_toggle(self):
        state = "normal" if self.mod_darkmatter.get() else "disabled"
        if hasattr(self, "sp_rdrp_minscore"):
            self.sp_rdrp_minscore.config(state=state)
        for btn in getattr(self, "dm_preset_buttons", []):
            btn.config(state=state)

    # =========================================================================
    # TAB 3: REAL-TIME TERMINAL CONSOLE BUILDER
    # =========================================================================
    def build_tab_terminal(self, parent):
        """
        Builds the contents of the Terminal Console tab with enlarged monospace font.
        """
        container = ttk.Frame(parent)
        container.pack(fill="both", expand=True, padx=6, pady=6)

        # Control bar on top of the console
        ctrl_bar = tk.Frame(container, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        ctrl_bar.pack(fill="x", pady=(0, 6))

        lbl_term_title = tk.Label(ctrl_bar, text="DeepVir Live Execution Console", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_term_title, "head", "bold")
        lbl_term_title.pack(side="left", padx=14, pady=10)

        # Live line counter
        self.lbl_lines = tk.Label(ctrl_bar, text="0 lines", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(self.lbl_lines, "body")
        self.lbl_lines.pack(side="left", padx=10)

        # Live stage & rule tracker badge
        self.stage_var = tk.StringVar(value="Stage: Idle")
        self.lbl_stage = tk.Label(ctrl_bar, textvariable=self.stage_var, fg=COLOR_ACCENT, bg="#0f172a", padx=12, pady=4)
        self.register_scalable(self.lbl_stage, "body", "bold")
        self.lbl_stage.pack(side="left", padx=8)

        # Auto-scroll toggle
        self.var_autoscroll = tk.BooleanVar(value=True)
        chk_scroll = ttk.Checkbutton(ctrl_bar, text="Auto-scroll", variable=self.var_autoscroll)
        chk_scroll.pack(side="right", padx=12)

        # Save log button
        btn_save_log = ttk.Button(ctrl_bar, text="[Save Log...]", style="Secondary.TButton", command=self.save_log_file)
        btn_save_log.pack(side="right", padx=6)

        # Clear console buffer button
        btn_clear = ttk.Button(ctrl_bar, text="[Clear]", style="Secondary.TButton", command=self.clear_terminal)
        btn_clear.pack(side="right", padx=6)

        # Scrolled Text Terminal with dark background and enlarged monospace text
        self.txt_terminal = scrolledtext.ScrolledText(
            container,
            bg=COLOR_CONSOLE_BG,
            fg=COLOR_CONSOLE_TEXT,
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )
        self.register_scalable(self.txt_terminal, "mono")
        self.txt_terminal.pack(fill="both", expand=True)

        # Configure color tags for ANSI-free structured log highlighting
        self.txt_terminal.tag_config("info", foreground="#38bdf8")      # Sky blue for [INFO]
        self.txt_terminal.tag_config("warn", foreground="#facc15")      # Yellow for [WARNING]
        self.txt_terminal.tag_config("error", foreground="#f87171")     # Red for [ERROR]
        self.txt_terminal.tag_config("success", foreground="#4ade80")   # Bright green for OK / Snakemake stages
        self.txt_terminal.tag_config("rule", foreground="#c084fc")      # Purple for Snakemake `rule <name>:`

    # =========================================================================
    # TAB 4: BUILT-IN RESULTS EXPLORER & SMART READER (.txt / .tsv / .html / .fasta)
    # =========================================================================
    def build_tab_results(self, parent):
        """
        Builds an interactive Results Explorer & Built-In Smart Reader in Tab 4:
        - Left pane: filterable file browser of generated outputs in `results/` and `logs/`.
        - Right pane: smart document reader that automatically adapts to `.txt`, `.tsv`,
          `.json`, `.fasta`, `.log`, and `.html` reports, plus WSL-compatible external
          browser/Explorer launchers.
        """
        self.results_filter_var = tk.StringVar(value="all")
        self.results_file_items = []  # List of Path objects currently in Listbox
        self.selected_result_path = None
        self._reader_images = []      # Keep references to embedded tk.PhotoImage charts
        self._cached_scan_entries = []  # In-memory cache: list of (Path, prefix_label, root_path, sz, mtime)
        self._cached_scan_out_dir = None
        self._html_image_cache = {}   # Cache decoded tk.PhotoImage by (len, head_snip) for instant re-render

        container = ttk.Frame(parent)
        container.pack(fill="both", expand=True, padx=6, pady=6)

        # Top Filter & Folder Bar
        top_bar = tk.Frame(container, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        top_bar.pack(fill="x", pady=(0, 6))

        top_row1 = tk.Frame(top_bar, bg=COLOR_CARD)
        top_row1.pack(fill="x", padx=12, pady=(8, 4))

        lbl_title = tk.Label(top_row1, text="Results Explorer & Smart File Reader", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_title, "head", "bold")
        lbl_title.pack(side="left", padx=(0, 12))

        self.lbl_results_count = tk.Label(top_row1, text="0 files found", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(self.lbl_results_count, "small", "bold")
        lbl_results_count_widget = self.lbl_results_count
        lbl_results_count_widget.pack(side="left")

        btn_open_res_dir = ttk.Button(
            top_row1,
            text="[Open results/ Folder]",
            style="Secondary.TButton",
            command=self.open_output_folder,
        )
        btn_open_res_dir.pack(side="right", padx=4)

        btn_open_logs_dir = ttk.Button(
            top_row1,
            text="[Open logs/ Folder]",
            style="Secondary.TButton",
            command=lambda: self.open_folder_in_os(LOGS_DIR),
        )
        btn_open_logs_dir.pack(side="right", padx=4)

        btn_refresh = ttk.Button(
            top_row1,
            text="[Refresh List]",
            style="Primary.TButton",
            command=lambda: self.refresh_results_list(rescan_disk=True),
        )
        btn_refresh.pack(side="right", padx=4)

        # Filter Category Buttons Row
        top_row2 = tk.Frame(top_bar, bg=COLOR_CARD)
        top_row2.pack(fill="x", padx=12, pady=(2, 8))

        lbl_flt = tk.Label(top_row2, text="Filter by Type:", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(lbl_flt, "small", "bold")
        lbl_flt.pack(side="left", padx=(0, 8))

        self.filter_buttons = {}
        for f_key, f_label in [
            ("all", "All Outputs"),
            ("html", "HTML Reports (.html)"),
            ("txt_stats", "Summaries & Stats (.txt / .tsv / .json)"),
            ("lineage", "Lineages & Hits (.tsv / .megan)"),
            ("fasta", "Contigs & ORFs (.fasta)"),
            ("logs", "Execution Logs (logs/)"),
        ]:
            btn_f = tk.Button(
                top_row2,
                text=f_label,
                bg=COLOR_PRIMARY if f_key == "all" else "#334155",
                fg="#ffffff" if f_key == "all" else COLOR_TEXT,
                relief="flat",
                padx=10,
                pady=3,
                cursor="hand2",
                activebackground=COLOR_PRIMARY_HOVER,
                activeforeground="#ffffff",
                command=lambda k=f_key: self.set_results_filter(k),
            )
            self.register_scalable(btn_f, "small", "bold")
            btn_f.pack(side="left", padx=3)
            self.filter_buttons[f_key] = btn_f

        # Window / Pane Resize Preset Controls (Right side of top_row2)
        pane_ctrl_box = tk.Frame(top_row2, bg=COLOR_CARD)
        pane_ctrl_box.pack(side="right", padx=(6, 0))

        lbl_pane_sz = tk.Label(pane_ctrl_box, text="Pane Size:", fg=COLOR_MUTED, bg=COLOR_CARD)
        self.register_scalable(lbl_pane_sz, "small", "bold")
        lbl_pane_sz.pack(side="left", padx=(0, 4))

        self.results_list_hidden = False
        for preset_label, preset_mode in [
            ("[Hide/Show List]", "toggle_list"),
            ("[<< Max Reader]", "max_reader"),
            ("[Standard]", "standard"),
            ("[50 / 50]", "equal"),
            ("[Wide Files >>]", "wide_files"),
        ]:
            btn_p = tk.Button(
                pane_ctrl_box,
                text=preset_label,
                bg="#1e293b",
                fg=COLOR_ACCENT if preset_mode != "toggle_list" else "#facc15",
                relief="flat",
                padx=7,
                pady=2,
                cursor="hand2",
                activebackground=COLOR_PRIMARY,
                activeforeground="#ffffff",
                command=lambda m=preset_mode: self.set_results_pane_layout(m),
            )
            self.register_scalable(btn_p, "small", "bold")
            btn_p.pack(side="left", padx=2)

        # Main Draggable Split Area: Left (File List) <--> Right (Smart Reader)
        self.results_paned = tk.PanedWindow(
            container,
            orient="horizontal",
            bg="#334155",
            sashwidth=8,
            sashpad=2,
            sashrelief="raised",
            opaqueresize=True,
            bd=0,
        )
        self.results_paned.pack(fill="both", expand=True)

        # Left Pane: Resizable File Browser List
        self.results_left_pane = tk.Frame(
            self.results_paned,
            bg=COLOR_CARD,
            width=430,
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )

        left_hdr_row = tk.Frame(self.results_left_pane, bg=COLOR_CARD)
        left_hdr_row.pack(fill="x", padx=10, pady=(8, 2))

        lbl_list_hdr = tk.Label(
            left_hdr_row,
            text="Output Files (drag divider <-> to resize):",
            fg=COLOR_PRIMARY,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_list_hdr, "body", "bold")
        lbl_list_hdr.pack(side="left")

        self.var_selected_relpath = tk.StringVar(value="Path: (select a file)")
        lbl_sel_path = tk.Label(
            self.results_left_pane,
            textvariable=self.var_selected_relpath,
            fg=COLOR_ACCENT,
            bg="#0f172a",
            anchor="w",
            padx=8,
            pady=3,
        )
        self.register_scalable(lbl_sel_path, "mono")
        lbl_sel_path.pack(fill="x", padx=8, pady=(0, 4))

        list_box_frame = tk.Frame(self.results_left_pane, bg=COLOR_CARD)
        list_box_frame.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        list_vscroll = ttk.Scrollbar(list_box_frame, orient="vertical")
        list_hscroll = ttk.Scrollbar(list_box_frame, orient="horizontal")

        self.lst_results = tk.Listbox(
            list_box_frame,
            bg="#0b0f19",
            fg=COLOR_TEXT,
            selectbackground=COLOR_PRIMARY,
            selectforeground="#ffffff",
            activestyle="none",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
            yscrollcommand=list_vscroll.set,
            xscrollcommand=list_hscroll.set,
        )
        self.register_scalable(self.lst_results, "mono")
        list_vscroll.config(command=self.lst_results.yview)
        list_hscroll.config(command=self.lst_results.xview)

        list_vscroll.pack(side="right", fill="y")
        list_hscroll.pack(side="bottom", fill="x")
        self.lst_results.pack(side="left", fill="both", expand=True)
        self.lst_results.bind("<<ListboxSelect>>", self.on_result_file_selected)
        self.lst_results.bind("<Double-Button-1>", lambda _e: self.open_selected_in_external_app())

        # Right Pane: Built-In Smart Reader
        self.results_right_pane = tk.Frame(
            self.results_paned,
            bg=COLOR_CARD,
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
        )

        self.results_paned.add(self.results_left_pane, minsize=180, width=430, stretch="never")
        self.results_paned.add(self.results_right_pane, minsize=320, stretch="always")
        right_pane = self.results_right_pane

        reader_bar = tk.Frame(right_pane, bg="#172033", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        reader_bar.pack(fill="x", padx=8, pady=(8, 4))

        self.var_reader_badge = tk.StringVar(value="[READER IDLE]")
        self.lbl_reader_badge = tk.Label(
            reader_bar,
            textvariable=self.var_reader_badge,
            fg="#10b981",
            bg="#064e3b",
            padx=10,
            pady=4,
        )
        self.register_scalable(self.lbl_reader_badge, "small", "bold")
        self.lbl_reader_badge.pack(side="left", padx=8, pady=6)

        self.var_reader_title = tk.StringVar(value="Select any file on the left to read its content here")
        lbl_reader_file = tk.Label(reader_bar, textvariable=self.var_reader_title, fg=COLOR_TEXT, bg="#172033")
        self.register_scalable(lbl_reader_file, "body", "bold")
        lbl_reader_file.pack(side="left", padx=6, pady=6)

        self.btn_open_external = ttk.Button(
            reader_bar,
            text="[Open in Browser / External App]",
            style="Primary.TButton",
            command=self.open_selected_in_external_app,
        )
        self.btn_open_external.pack(side="right", padx=8, pady=6)

        # Reader Font Zoom & Line Wrap Controls
        btn_font_plus = tk.Button(
            reader_bar,
            text="[Font +]",
            bg="#334155",
            fg="#ffffff",
            relief="flat",
            padx=7,
            pady=2,
            cursor="hand2",
            command=lambda: self.adjust_scale(1),
        )
        self.register_scalable(btn_font_plus, "small", "bold")
        btn_font_plus.pack(side="right", padx=2, pady=6)

        btn_font_minus = tk.Button(
            reader_bar,
            text="[Font -]",
            bg="#334155",
            fg="#ffffff",
            relief="flat",
            padx=7,
            pady=2,
            cursor="hand2",
            command=lambda: self.adjust_scale(-1),
        )
        self.register_scalable(btn_font_minus, "small", "bold")
        btn_font_minus.pack(side="right", padx=2, pady=6)

        self.var_reader_wrap_btn = tk.StringVar(value="[Wrap: AUTO]")
        btn_toggle_wrap = tk.Button(
            reader_bar,
            textvariable=self.var_reader_wrap_btn,
            bg="#334155",
            fg=COLOR_ACCENT,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.toggle_reader_line_wrap,
        )
        self.register_scalable(btn_toggle_wrap, "small", "bold")
        btn_toggle_wrap.pack(side="right", padx=4, pady=6)

        # Dark Matter (RdRp & RVDB) Hit-to-FASTA Linker & Threshold Extractor Bar
        dm_bar = tk.Frame(right_pane, bg="#0f172a", highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        dm_bar.pack(fill="x", padx=8, pady=(0, 4))

        lbl_dm_bar = tk.Label(
            dm_bar,
            text="Dark Matter Linker:",
            fg=COLOR_ACCENT,
            bg="#0f172a",
        )
        self.register_scalable(lbl_dm_bar, "small", "bold")
        lbl_dm_bar.pack(side="left", padx=(8, 4), pady=5)

        lbl_dm_ms = tk.Label(dm_bar, text="RdRp MinScore >=", fg=COLOR_TEXT, bg="#0f172a")
        self.register_scalable(lbl_dm_ms, "small")
        lbl_dm_ms.pack(side="left", padx=(4, 2), pady=5)

        init_ms = 75
        if hasattr(self, "var_rdrp_minscore"):
            try:
                init_ms = int(self.var_rdrp_minscore.get())
            except (ValueError, tk.TclError):
                init_ms = 75
        self.var_dm_filter_minscore = tk.IntVar(value=init_ms)
        self.sp_dm_filter_minscore = tk.Spinbox(
            dm_bar,
            from_=0,
            to=100,
            increment=5,
            textvariable=self.var_dm_filter_minscore,
            width=5,
            bg="#1e293b",
            fg=COLOR_TEXT,
            insertbackground="white",
        )
        self.register_scalable(self.sp_dm_filter_minscore, "mono", "bold")
        self.sp_dm_filter_minscore.pack(side="left", padx=(0, 8), pady=5)

        lbl_rvdb_ev = tk.Label(dm_bar, text="RVDB E-val <=", fg=COLOR_TEXT, bg="#0f172a")
        self.register_scalable(lbl_rvdb_ev, "small")
        lbl_rvdb_ev.pack(side="left", padx=(4, 2), pady=5)

        self.var_rvdb_max_evalue = tk.StringVar(value="1e-5")
        self.cb_rvdb_evalue = ttk.Combobox(
            dm_bar,
            textvariable=self.var_rvdb_max_evalue,
            values=["1e-20", "1e-10", "1e-5", "1e-3", "1.0"],
            width=7,
        )
        self.register_scalable(self.cb_rvdb_evalue, "small")
        self.cb_rvdb_evalue.pack(side="left", padx=(0, 8), pady=5)

        self.var_rvdb_cat = tk.StringVar(value="All (Pol + Struct)")
        self.cb_rvdb_cat = ttk.Combobox(
            dm_bar,
            textvariable=self.var_rvdb_cat,
            values=["All (Pol + Struct)", "Polymerase / RdRp", "Structural Only"],
            width=16,
            state="readonly",
        )
        self.register_scalable(self.cb_rvdb_cat, "small")
        self.cb_rvdb_cat.pack(side="left", padx=(0, 8), pady=5)

        btn_extract_fastas = tk.Button(
            dm_bar,
            text="[Extract & Save Linked FASTAs]",
            bg="#059669",
            fg="#ffffff",
            relief="flat",
            padx=10,
            pady=3,
            cursor="hand2",
            activebackground="#10b981",
            activeforeground="#ffffff",
            command=self.extract_and_save_linked_fastas,
        )
        self.register_scalable(btn_extract_fastas, "small", "bold")
        btn_extract_fastas.pack(side="right", padx=(4, 8), pady=4)

        btn_preview_linked = tk.Button(
            dm_bar,
            text="[Link & View Hits >= Threshold]",
            bg=COLOR_PRIMARY,
            fg="#ffffff",
            relief="flat",
            padx=10,
            pady=3,
            cursor="hand2",
            activebackground=COLOR_PRIMARY_HOVER,
            activeforeground="#ffffff",
            command=self.preview_linked_darkmatter_hits,
        )
        self.register_scalable(btn_preview_linked, "small", "bold")
        btn_preview_linked.pack(side="right", padx=4, pady=4)

        # Smart Reader ScrolledText
        self.txt_reader = scrolledtext.ScrolledText(
            right_pane,
            bg=COLOR_CONSOLE_BG,
            fg=COLOR_TEXT,
            insertbackground="white",
            relief="flat",
            highlightbackground=COLOR_CARD_BORDER,
            highlightthickness=1,
            wrap="none",
        )
        self.register_scalable(self.txt_reader, "mono")
        self.txt_reader.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        reader_hscroll = ttk.Scrollbar(right_pane, orient="horizontal", command=self.txt_reader.xview)
        self.txt_reader.configure(xscrollcommand=reader_hscroll.set)
        reader_hscroll.pack(fill="x", padx=8, pady=(0, 6))

        # Configure color tags for the Smart Reader
        self.txt_reader.tag_config("h1", foreground=COLOR_PRIMARY)
        self.txt_reader.tag_config("h2", foreground=COLOR_ACCENT)
        self.txt_reader.tag_config("h3", foreground="#a78bfa")
        self.txt_reader.tag_config("bold", foreground="#ffffff")
        self.txt_reader.tag_config("table_hdr", foreground="#facc15")
        self.txt_reader.tag_config("table_border", foreground="#475569")
        self.txt_reader.tag_config("bar_fill", foreground="#10b981")
        self.txt_reader.tag_config("fasta_hdr", foreground="#38bdf8")
        self.txt_reader.tag_config("key", foreground="#34d399")
        self.txt_reader.tag_config("code", foreground="#fde047", background="#1e293b")
        self.txt_reader.tag_config("muted", foreground=COLOR_MUTED)
        self.txt_reader.tag_config("error", foreground="#f87171")

        # Initial scan
        self.after(200, lambda: self.refresh_results_list(rescan_disk=True))

    def _on_notebook_tab_changed(self, _event=None):
        """Refreshes the Results Explorer when switching to Tab 4 (uses fast cache if already scanned)."""
        try:
            if self.notebook.index(self.notebook.select()) == 3:
                out_dir_str = self.var_output.get().strip() or str(PROJECT_DIR / "results")
                need_rescan = (not self._cached_scan_entries) or (self._cached_scan_out_dir != out_dir_str)
                self.refresh_results_list(rescan_disk=need_rescan)
        except Exception:
            pass

    def set_results_pane_layout(self, mode: str):
        """
        Resizes or toggles the left file browser and right reader panes in Tab 4:
        - 'toggle_list': hides/restores the left file list so the reader uses 100% width
        - 'max_reader': narrows left pane to 220px to maximize reader space
        - 'standard': sets left pane to 430px
        - 'equal': splits 50% / 50% between file list and reader
        - 'wide_files': expands left pane to 62% width for reading long file paths
        """
        if not hasattr(self, "results_paned"):
            return
        self.update_idletasks()
        total_w = max(700, self.results_paned.winfo_width())

        if mode == "toggle_list":
            panes = self.results_paned.panes()
            if str(self.results_left_pane) in panes:
                self.results_paned.forget(self.results_left_pane)
                self.results_list_hidden = True
            else:
                self.results_paned.add(
                    self.results_left_pane,
                    before=self.results_right_pane,
                    minsize=180,
                    width=min(430, int(total_w * 0.35)),
                    stretch="never",
                )
                self.results_list_hidden = False
            return

        # Ensure left pane is visible before moving sash
        if str(self.results_left_pane) not in self.results_paned.panes():
            self.results_paned.add(
                self.results_left_pane,
                before=self.results_right_pane,
                minsize=180,
                width=430,
                stretch="never",
            )
            self.results_list_hidden = False
            self.update_idletasks()

        if mode == "max_reader":
            target_x = 220
        elif mode == "standard":
            target_x = min(430, int(total_w * 0.36))
        elif mode == "equal":
            target_x = int(total_w * 0.50)
        elif mode == "wide_files":
            target_x = int(total_w * 0.62)
        else:
            target_x = 430

        try:
            self.results_paned.sash_place(0, target_x, 1)
        except Exception:
            pass

    def toggle_reader_line_wrap(self):
        """Toggles horizontal line wrapping (`none` <-> `word`) inside `self.txt_reader`."""
        if not hasattr(self, "txt_reader"):
            return
        cur_wrap = str(self.txt_reader.cget("wrap"))
        if cur_wrap == "none":
            self.txt_reader.configure(wrap="word")
            if hasattr(self, "var_reader_wrap_btn"):
                self.var_reader_wrap_btn.set("[Wrap: ON]")
        else:
            self.txt_reader.configure(wrap="none")
            if hasattr(self, "var_reader_wrap_btn"):
                self.var_reader_wrap_btn.set("[Wrap: OFF]")

    def set_results_filter(self, filter_key: str):
        self.results_filter_var.set(filter_key)
        for k, btn in self.filter_buttons.items():
            if k == filter_key:
                btn.config(bg=COLOR_PRIMARY, fg="#ffffff")
            else:
                btn.config(bg="#334155", fg=COLOR_TEXT)
        # Filter in-memory cache instantaneously without touching disk
        self.refresh_results_list(rescan_disk=False)

    def _scan_results_disk(self, out_dir: Path) -> list:
        """
        Fast directory scanner using `os.walk` with aggressive directory pruning so
        heavy assembler intermediate folders (`K21`, `K33`, `tmp`, `intermediate_contigs`,
        `multiqc_data`, `.snakemake`) are skipped immediately without stat() overhead.
        """
        READABLE_EXTS = {
            ".html", ".htm", ".txt", ".tsv", ".tab", ".csv",
            ".json", ".fasta", ".fa", ".fas", ".fna", ".megan",
            ".log", ".yaml", ".yml", ".fev",
        }
        SKIP_DIRS = {
            "benchmarks", ".snakemake", "intermediate_contigs", "multiqc_data",
            "multiqc_report_data", "all_samples_multiqc_aggregate_data",
            "K21", "K33", "K55", "K77", "K99", "K127",
            "corrected", "tmp", "misc", "split_input", "__pycache__",
            "configs", "pipeline_state",
        }

        entries = []
        for root_path, prefix_label in [(out_dir, "results"), (LOGS_DIR, "logs")]:
            if not root_path.is_dir():
                continue
            for dirpath, dirnames, filenames in os.walk(root_path):
                # Prune heavy intermediate directories in-place before descending
                dirnames[:] = [
                    d for d in dirnames
                    if d not in SKIP_DIRS
                    and not d.startswith(".")
                    and not d.endswith("_multiqc_sample_data")
                ]
                for fname in filenames:
                    ext = os.path.splitext(fname)[1].lower()
                    if ext not in READABLE_EXTS:
                        continue
                    p = Path(dirpath) / fname
                    try:
                        st = p.stat()
                        sz = st.st_size
                        mtime = st.st_mtime
                    except OSError:
                        sz = 0
                        mtime = 0.0
                    entries.append((p, prefix_label, root_path, sz, mtime))

        # Sort once: HTML & summary reports first, then by modification time (newest first)
        def _sort_key(item):
            p, prefix_label, _, _, mtime = item
            ext = p.suffix.lower()
            priority = 0 if (prefix_label == "results" and ext in (".html", ".txt", ".tsv", ".csv", ".json")) else 1
            return (priority, -mtime, str(p))

        entries.sort(key=_sort_key)
        return entries

    def refresh_results_list(self, rescan_disk: bool = True):
        """Populates Tab 4's file list from the fast pruned scanner or in-memory cache."""
        if not hasattr(self, "lst_results"):
            return

        flt = self.results_filter_var.get()
        out_dir_str = self.var_output.get().strip() or str(PROJECT_DIR / "results")
        out_dir = Path(out_dir_str)

        if rescan_disk or not self._cached_scan_entries or self._cached_scan_out_dir != out_dir_str:
            self._cached_scan_entries = self._scan_results_disk(out_dir)
            self._cached_scan_out_dir = out_dir_str

        candidates = []
        for p, prefix_label, root_path, sz, mtime in self._cached_scan_entries:
            if flt == "logs" and prefix_label != "logs":
                continue
            if flt not in ("all", "logs") and prefix_label != "results":
                continue

            ext = p.suffix.lower()
            name_low = p.name.lower()

            if flt == "html" and ext not in (".html", ".htm"):
                continue
            elif flt == "txt_stats" and not (
                ext in (".txt", ".json")
                or "stats" in name_low
                or "summary" in name_low
                or "report" in name_low
                or "biom" in name_low
                or "diamond_view" in name_low
                or "kraken2_" in name_low
            ):
                continue
            elif flt == "lineage" and not (
                "lineage" in name_low
                or "hits" in name_low
                or "rdrp" in name_low
                or "fev2tsv" in name_low
                or "palm_annot" in name_low
                or "rvdb" in name_low
                or "daa2info" in name_low
                or ext == ".megan"
                or "otu_table" in name_low
            ):
                continue
            elif flt == "fasta" and ext not in (".fasta", ".fa", ".fas", ".fna"):
                continue
            elif flt == "logs" and ext not in (".log", ".yaml", ".yml", ".txt"):
                continue

            candidates.append((p, prefix_label, root_path, sz))

        self.lst_results.delete(0, "end")
        self.results_file_items = []
        display_lines = []

        for p, prefix_label, root_path, sz in candidates:
            try:
                rel = p.relative_to(root_path)
            except ValueError:
                rel = Path(p.name)
            if sz < 1024:
                sz_str = f"{sz} B"
            elif sz < 1024 * 1024:
                sz_str = f"{sz / 1024:.1f} KB"
            else:
                sz_str = f"{sz / (1024 * 1024):.1f} MB"

            ext_tag = p.suffix.upper().lstrip(".") or "FILE"
            display_lines.append(f"[{ext_tag:<5}] {prefix_label}/{rel}  ({sz_str})")
            self.results_file_items.append(p)

        if display_lines:
            self.lst_results.insert("end", *display_lines)

        self.lbl_results_count.config(text=f"{len(self.results_file_items)} files found")

        # Auto-select first file if nothing is currently loaded
        if self.results_file_items and (
            self.selected_result_path not in self.results_file_items
        ):
            self.lst_results.selection_set(0)
            self.load_file_into_reader(self.results_file_items[0])
        elif not self.results_file_items:
            self.txt_reader.delete("1.0", "end")
            self.txt_reader.insert(
                "end",
                "No output files found for this filter yet.\n\n"
                "Run the pipeline (or change the filter above) and click [Refresh List].\n",
                "muted",
            )

    def on_result_file_selected(self, _event=None):
        sel = self.lst_results.curselection()
        if not sel:
            return
        idx = sel[0]
        if 0 <= idx < len(self.results_file_items):
            self.load_file_into_reader(self.results_file_items[idx])

    def load_file_into_reader(self, file_path: Path):
        """
        Smart In-App Reader:
        - If `file_path` is `.html` / `.htm`, renders a spacious visual document directly
          inside `self.txt_reader` with word wrapping, boxed tables, visual bar charts
          (MultiQC / Krona), and embedded PNG figures (`tk.PhotoImage`).
        - If `file_path` is `.tsv`/`.tab`, `.json`, `.fasta`/`.fa`, or `.txt`/`.log`,
          formats it for clean reading inside the built-in GUI reader.
        """
        self.selected_result_path = file_path
        self._reader_images.clear()
        ext = file_path.suffix.lower()

        try:
            sz = file_path.stat().st_size
        except OSError:
            sz = 0

        self.var_reader_title.set(f"{file_path.name}  ({sz:,} bytes)")
        if hasattr(self, "var_selected_relpath"):
            out_dir = Path(self.var_output.get().strip() or str(PROJECT_DIR / "results"))
            try:
                rel_disp = f"results/{file_path.relative_to(out_dir)}"
            except ValueError:
                try:
                    rel_disp = f"logs/{file_path.relative_to(LOGS_DIR)}"
                except ValueError:
                    rel_disp = str(file_path.name)
            self.var_selected_relpath.set(f"Selected: {rel_disp}  ({sz:,} B)")
        self.txt_reader.delete("1.0", "end")

        if ext in (".html", ".htm"):
            # Spacious document layout for HTML reports inside the app
            self.txt_reader.configure(wrap="word", padx=20, pady=14, spacing1=4, spacing2=2, spacing3=5)
            if hasattr(self, "var_reader_wrap_btn"):
                self.var_reader_wrap_btn.set("[Wrap: ON]")
            self.var_reader_badge.set("[HTML VISUAL REPORT]")
            MAX_HTML_BYTES = 30 * 1024 * 1024  # Read up to 30 MB so embedded base64 PNG plots load completely
            try:
                with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                    raw_html = f.read(MAX_HTML_BYTES)
            except Exception as e:
                self.var_reader_badge.set("[READ ERROR]")
                self.txt_reader.insert("end", f"Could not read HTML file {file_path}:\n{e}\n", "error")
                return
            self._render_html_in_reader(file_path, raw_html)
            return

        # Standard monospace layout for tables, FASTA, JSON, and logs
        self.txt_reader.configure(wrap="none", padx=12, pady=10, spacing1=1, spacing2=0, spacing3=1)
        if hasattr(self, "var_reader_wrap_btn"):
            self.var_reader_wrap_btn.set("[Wrap: OFF]")

        # Read up to 128 KB for instant preview of large text/table/FASTA files
        MAX_BYTES = 128 * 1024
        try:
            with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                raw_content = f.read(MAX_BYTES)
            truncated = sz > MAX_BYTES
        except Exception as e:
            self.var_reader_badge.set("[READ ERROR]")
            self.txt_reader.insert("end", f"Could not read file {file_path}:\n{e}\n", "error")
            return

        if ext in (".tsv", ".tab") or (ext == ".csv" and "\t" in raw_content[:2048]):
            self.var_reader_badge.set("[TSV TABLE]")
            self._render_tsv_in_reader(file_path, raw_content, truncated)
        elif ext == ".json":
            self.var_reader_badge.set("[JSON DATA]")
            self._render_text_in_reader(file_path, raw_content, truncated)
        elif ext in (".fasta", ".fa", ".fas", ".fna"):
            self.var_reader_badge.set("[FASTA SEQUENCES]")
            self._render_fasta_in_reader(file_path, raw_content, truncated)
        else:
            self.var_reader_badge.set("[TEXT / LOG]")
            self._render_text_in_reader(file_path, raw_content, truncated)

    def _render_text_in_reader(self, file_path: Path, content: str, truncated: bool):
        self.txt_reader.insert("end", f"File: {file_path}\n", "muted")
        self.txt_reader.insert("end", "=" * 78 + "\n\n", "h2")

        for line in content.splitlines(keepends=True):
            stripped = line.strip()
            if stripped.startswith("===") or stripped.startswith("---"):
                self.txt_reader.insert("end", line, "h2")
            elif stripped.startswith("[ERROR]") or "Error" in stripped:
                self.txt_reader.insert("end", line, "error")
            elif ":" in stripped and not stripped.startswith("http"):
                k, v = line.split(":", 1)
                self.txt_reader.insert("end", k + ":", "key")
                self.txt_reader.insert("end", v)
            else:
                self.txt_reader.insert("end", line)

        if truncated:
            self.txt_reader.insert(
                "end",
                "\n\n[NOTE: File preview truncated at 256 KB. Click 'Open in Browser / External App' to view full file.]\n",
                "muted",
            )

    def _render_tsv_in_reader(self, file_path: Path, content: str, truncated: bool):
        self.txt_reader.insert("end", f"Table File: {file_path}\n", "muted")
        lines = [ln for ln in content.splitlines() if ln.strip()]
        if not lines:
            self.txt_reader.insert("end", "(Empty table file - 0 rows)\n", "muted")
            return

        rows = [ln.split("\t") for ln in lines[:300]]
        header = [c.strip() for c in rows[0]] if rows else []

        # Detect Dark Matter RdRp TSV (`_RdRp.tsv`) or RVDB Summary (`_RVDB_Summary.csv` / `_RVDB_results.tsv`)
        is_rdrp_tsv = "Label" in header and "rdrp" in header
        is_rvdb_tsv = "Sequence_ID" in header and "E_value" in header

        if is_rdrp_tsv or is_rvdb_tsv:
            try:
                ms_thresh = float(self.var_dm_filter_minscore.get())
            except (ValueError, tk.TclError):
                ms_thresh = 75.0
            try:
                ev_thresh = float(self.var_rvdb_max_evalue.get())
            except (ValueError, tk.TclError):
                ev_thresh = 1e-5

            # Prepend Linked_Contig_ID and Pass_Threshold columns dynamically for instant linking
            label_idx = header.index("Label") if is_rdrp_tsv else header.index("Sequence_ID")
            score_idx = header.index("rdrp") if is_rdrp_tsv else header.index("E_value")

            augmented_rows = [["Status", "Linked_Contig_ID"] + header]
            pass_count = 0
            linked_contigs = set()

            for r in rows[1:]:
                orf_label = r[label_idx].strip() if label_idx < len(r) else ""
                contig_id = self._extract_contig_id_from_orf_label(orf_label)
                raw_score = r[score_idx].strip() if score_idx < len(r) else ""
                passes = False
                try:
                    val_f = float(raw_score)
                    passes = (val_f >= ms_thresh) if is_rdrp_tsv else (val_f <= ev_thresh)
                except ValueError:
                    passes = False
                if passes:
                    pass_count += 1
                    if contig_id != "-":
                        linked_contigs.add(contig_id)
                status_tag = "[>= THRESH]" if passes else "[below]"
                augmented_rows.append([status_tag, contig_id] + r)

            rows = augmented_rows
            if is_rdrp_tsv:
                self.var_reader_badge.set("[RdRp DARK MATTER HITS]")
                self.txt_reader.insert(
                    "end",
                    f"Dark Matter RdRp Hits: {pass_count}/{max(0, len(lines) - 1)} ORFs >= MinScore ({ms_thresh:g}) "
                    f"| Linked to {len(linked_contigs)} unique Contigs\n",
                    "bar_fill",
                )
            else:
                self.var_reader_badge.set("[RVDB VALIDATION HITS]")
                self.txt_reader.insert(
                    "end",
                    f"RVDB Validation Hits: {pass_count}/{max(0, len(lines) - 1)} ORFs <= E-value ({ev_thresh:g}) "
                    f"| Linked to {len(linked_contigs)} unique Contigs\n",
                    "bar_fill",
                )
            self.txt_reader.insert(
                "end",
                "TIP: Adjust threshold in the 'Dark Matter Linker' bar above and click "
                "[Link & View Hits >= Threshold] or [Extract & Save Linked FASTAs]!\n",
                "h3",
            )

        num_cols = max(len(r) for r in rows)
        col_widths = [0] * num_cols
        for r in rows:
            for c_idx, cell in enumerate(r):
                col_widths[c_idx] = min(42, max(col_widths[c_idx], len(cell.strip())))

        total_rows = max(0, len(lines) - 1)
        self.txt_reader.insert("end", f"Rows: {total_rows} data rows | Columns: {num_cols}\n", "h1")
        self.txt_reader.insert("end", "=" * 88 + "\n", "h2")

        # Batch consecutive rows sharing the same tag to minimize Tkinter Tcl round-trips
        cur_tag = None
        cur_buf = []
        for r_idx, r in enumerate(rows[:150]):
            formatted_cells = []
            for c_idx in range(num_cols):
                val = r[c_idx].strip() if c_idx < len(r) else ""
                if len(val) > 42:
                    val = val[:39] + "..."
                formatted_cells.append(val.ljust(col_widths[c_idx]))
            row_str = " | ".join(formatted_cells) + "\n"
            if r_idx == 0:
                self.txt_reader.insert("end", row_str, "table_hdr")
                sep = "-+-".join("-" * w for w in col_widths) + "\n"
                self.txt_reader.insert("end", sep, "h2")
            else:
                if (is_rdrp_tsv or is_rvdb_tsv) and r and r[0] == "[>= THRESH]":
                    row_tag = "bar_fill"
                elif (is_rdrp_tsv or is_rvdb_tsv) and r and r[0] == "[below]":
                    row_tag = "muted"
                else:
                    row_tag = None
                if row_tag != cur_tag and cur_buf:
                    if cur_tag:
                        self.txt_reader.insert("end", "".join(cur_buf), cur_tag)
                    else:
                        self.txt_reader.insert("end", "".join(cur_buf))
                    cur_buf.clear()
                cur_tag = row_tag
                cur_buf.append(row_str)

        if cur_buf:
            if cur_tag:
                self.txt_reader.insert("end", "".join(cur_buf), cur_tag)
            else:
                self.txt_reader.insert("end", "".join(cur_buf))

        if len(lines) > 150 or truncated:
            self.txt_reader.insert(
                "end",
                "\n[Showing first 150 rows. Click 'Open in Browser / External App' for full table.]\n",
                "muted",
            )

    def _render_fasta_in_reader(self, file_path: Path, content: str, truncated: bool):
        seq_count = 0
        total_bp = 0
        seq_buf = []
        for ln in content.splitlines(keepends=True):
            if ln.startswith(">"):
                seq_count += 1
            else:
                total_bp += len(ln.strip())

        self.txt_reader.insert("end", f"FASTA File: {file_path}\n", "muted")
        self.txt_reader.insert(
            "end",
            f"Preview Stats: {seq_count} sequences shown | {total_bp:,} bp in preview\n",
            "h1",
        )
        self.txt_reader.insert("end", "=" * 78 + "\n\n", "h2")

        for ln in content.splitlines(keepends=True):
            if ln.startswith(">"):
                if seq_buf:
                    self.txt_reader.insert("end", "".join(seq_buf))
                    seq_buf.clear()
                self.txt_reader.insert("end", ln, "fasta_hdr")
            else:
                seq_buf.append(ln)
        if seq_buf:
            self.txt_reader.insert("end", "".join(seq_buf))

        if truncated:
            self.txt_reader.insert(
                "end",
                "\n\n[NOTE: FASTA preview truncated at 128 KB.]\n",
                "muted",
            )

    def _insert_formatted_box_table(self, caption: str, rows: list):
        """Renders a 2D table of strings as a clean, padded ASCII box table inside `self.txt_reader`."""
        if not rows:
            return
        num_cols = max(len(r) for r in rows)
        if num_cols == 0:
            return

        col_widths = [0] * num_cols
        for r in rows:
            for c_idx in range(num_cols):
                val = str(r[c_idx]).strip() if c_idx < len(r) else ""
                col_widths[c_idx] = min(44, max(col_widths[c_idx], len(val)))

        if caption:
            self.txt_reader.insert("end", f"\n  TABLE: {caption}\n", "h2")
        else:
            self.txt_reader.insert("end", "\n")

        border_line = "  +" + "+".join("-" * (w + 2) for w in col_widths) + "+\n"
        self.txt_reader.insert("end", border_line, "table_border")

        body_lines = []
        for r_idx, r in enumerate(rows):
            cells = []
            for c_idx in range(num_cols):
                val = str(r[c_idx]).strip() if c_idx < len(r) else ""
                if len(val) > 44:
                    val = val[:41] + "..."
                is_num = r_idx > 0 and bool(re.match(r"^-?[\d,]+(\.\d+)?%?$", val))
                padded = val.rjust(col_widths[c_idx]) if is_num else val.ljust(col_widths[c_idx])
                cells.append(f" {padded} ")

            row_line = "  |" + "|".join(cells) + "|\n"
            if r_idx == 0:
                self.txt_reader.insert("end", row_line, "table_hdr")
                self.txt_reader.insert("end", border_line, "table_border")
            else:
                body_lines.append(row_line)

        if body_lines:
            self.txt_reader.insert("end", "".join(body_lines), "bold")

        self.txt_reader.insert("end", border_line, "table_border")
        self.txt_reader.insert("end", "\n")

    def _insert_base64_png_image(self, src_data: str):
        """Decodes a `data:image/png;base64,...` URI (with in-memory caching) and embeds it in `self.txt_reader`."""
        try:
            b64_part = src_data.split("base64,", 1)[1].strip()
            cache_key = (len(b64_part), b64_part[:96], b64_part[-96:])
            photo = self._html_image_cache.get(cache_key)
            if photo is None:
                photo = tk.PhotoImage(data=b64_part)
                w = photo.width()
                target_w = max(580, self.txt_reader.winfo_width() - 80) if self.txt_reader.winfo_width() > 300 else 680
                if w > target_w:
                    scale = max(2, int(round(w / target_w)))
                    photo = photo.subsample(scale, scale)
                if len(self._html_image_cache) >= 32:
                    self._html_image_cache.clear()
                self._html_image_cache[cache_key] = photo
            self._reader_images.append(photo)
            self.txt_reader.insert("end", "\n  ")
            self.txt_reader.image_create("end", image=photo, padx=10, pady=8)
            self.txt_reader.insert("end", "\n\n")
        except Exception:
            pass

    def _render_multiqc_compressed_data(self, content: str) -> bool:
        """
        Extracts and decompresses MultiQC's `<script id="mqc_compressed_plotdata">`
        (base64 + gzip JSON) using fast O(1) string indexing instead of regex backtracking.
        """
        idx_id = content.find("mqc_compressed_plotdata")
        if idx_id == -1:
            return False
        idx_gt = content.find(">", idx_id)
        if idx_gt == -1:
            return False
        idx_end = content.find("</script>", idx_gt)
        if idx_end == -1:
            return False

        b64_str = content[idx_gt + 1 : idx_end].strip()
        try:
            raw_bytes = base64.b64decode(b64_str)
            json_str = gzip.decompress(raw_bytes).decode("utf-8", errors="replace")
            mqc = json.loads(json_str)
        except Exception:
            return False

        title = mqc.get("config_title") or "MultiQC Consolidated Quality Report"
        subtitle = mqc.get("config_subtitle") or ""
        created = mqc.get("config_creation_date") or ""

        self.txt_reader.insert("end", "=" * 84 + "\n", "h2")
        self.txt_reader.insert("end", f"  {title}\n", "h1")
        if subtitle:
            self.txt_reader.insert("end", f"  {subtitle}\n", "h3")
        if created:
            self.txt_reader.insert("end", f"  Generated: {created}\n", "muted")
        self.txt_reader.insert("end", "=" * 84 + "\n\n", "h2")

        # 1. Render General Statistics Table
        gen_stats = mqc.get("report_general_stats_data", [])
        gen_headers = mqc.get("report_general_stats_headers", [])
        if gen_stats:
            sample_metrics = {}
            metric_titles = []
            # MultiQC stores general stats either as list of dicts or dict of dicts
            stats_iter = gen_stats.values() if isinstance(gen_stats, dict) else gen_stats
            headers_iter = gen_headers.values() if isinstance(gen_headers, dict) else gen_headers

            hdr_lookup = {}
            for h_group in headers_iter:
                if isinstance(h_group, dict):
                    for k, h_meta in h_group.items():
                        if isinstance(h_meta, dict):
                            hdr_lookup[k] = h_meta.get("title") or k

            for s_group in stats_iter:
                if not isinstance(s_group, dict):
                    continue
                for sample_name, metrics_dict in s_group.items():
                    if not isinstance(metrics_dict, dict):
                        continue
                    sample_metrics.setdefault(sample_name, {})
                    for m_key, m_val in metrics_dict.items():
                        if isinstance(m_val, dict):
                            m_val = m_val.get("val", "")
                        if isinstance(m_val, (int, float)):
                            lbl = hdr_lookup.get(m_key, m_key)
                            if lbl not in metric_titles:
                                metric_titles.append(lbl)
                            if isinstance(m_val, float):
                                sample_metrics[sample_name][lbl] = f"{m_val:,.2f}"
                            else:
                                sample_metrics[sample_name][lbl] = f"{m_val:,}"

            if sample_metrics and metric_titles:
                shown_metrics = metric_titles[:8]
                table_rows = [["Sample"] + shown_metrics]
                for s_name, s_vals in sorted(sample_metrics.items()):
                    table_rows.append([s_name] + [s_vals.get(m_lbl, "-") for m_lbl in shown_metrics])
                self._insert_formatted_box_table("MultiQC General Statistics Summary", table_rows)

        # 2. Render Visual Bar Graphs & Plot Data from `report_plot_data`
        plot_data = mqc.get("report_plot_data", {})
        if isinstance(plot_data, dict):
            for plot_id, p_info in plot_data.items():
                if not isinstance(p_info, dict):
                    continue
                p_type = str(p_info.get("plot_type", "")).lower()
                p_cfg = p_info.get("config", {}) if isinstance(p_info.get("config"), dict) else {}
                p_title = p_cfg.get("title") or plot_id.replace("_", " ").title()
                datasets = p_info.get("datasets", [])

                if "bar" in p_type and datasets and isinstance(datasets[0], dict):
                    cats = datasets[0].get("cats", [])
                    samples = datasets[0].get("samples", [])
                    if cats and samples:
                        self.txt_reader.insert("end", f"  CHART: {p_title}\n", "h1")
                        self.txt_reader.insert("end", "  " + "-" * 78 + "\n", "h2")
                        for s_idx, s_name in enumerate(samples[:20]):
                            self.txt_reader.insert("end", f"  Sample: {s_name}\n", "table_hdr")
                            cat_vals = []
                            for c_dict in cats:
                                if not isinstance(c_dict, dict):
                                    continue
                                c_name = c_dict.get("name", "Category")
                                data_arr = c_dict.get("data", [])
                                val = data_arr[s_idx] if s_idx < len(data_arr) else 0
                                if isinstance(val, list) and val:
                                    val = val[0]
                                try:
                                    val_f = float(val)
                                except (TypeError, ValueError):
                                    val_f = 0.0
                                cat_vals.append((c_name, val_f))

                            tot = sum(v for _, v in cat_vals if v > 0) or 1.0
                            for c_name, val_f in cat_vals:
                                pct = (val_f / tot) * 100.0 if tot > 0 else 0.0
                                filled = int(round((pct / 100.0) * 24))
                                filled = max(0, min(24, filled))
                                bar_str = "[" + "#" * filled + "." * (24 - filled) + "]"
                                val_disp = f"{val_f:,.0f}" if val_f >= 10 else f"{val_f:,.2f}"
                                self.txt_reader.insert("end", f"    {c_name:<34} ")
                                self.txt_reader.insert("end", bar_str, "bar_fill")
                                self.txt_reader.insert("end", f"  {pct:5.1f}%  ({val_disp})\n", "muted")
                            self.txt_reader.insert("end", "\n")
        return True

    def _render_krona_hierarchy(self, content: str) -> bool:
        """Renders Krona `<node name="...">` XML taxonomic trees as a visual bar table in-app."""
        krona_idx = content.find("<krona")
        if krona_idx < 0:
            krona_idx = content.lower().find("<krona")
        if krona_idx < 0:
            return False

        krona_xml = content[krona_idx:]
        attr_m = re.search(r'<attributes[^>]*\bmagnitude="([^"]+)"', krona_xml, re.IGNORECASE)
        mag_tag = re.escape(attr_m.group(1)) if attr_m else "magnitude"
        val_re = re.compile(
            rf'<(?:{mag_tag}|count|magnitude)>\s*<val>\s*([+-]?\d+(?:\.\d+)?)\s*</val>',
            re.IGNORECASE,
        )
        rank_re = re.compile(r'<rank>\s*<val>\s*([^<]+?)\s*</val>', re.IGNORECASE)
        taxon_re = re.compile(r'<taxon>\s*<val[^>]*>\s*([^<]+?)\s*</val>', re.IGNORECASE)
        token_re = re.compile(r'(<node\s+name="([^"]+)"[^>]*>)|(</node\s*>)', re.IGNORECASE)

        matches = list(token_re.finditer(krona_xml))
        if not matches:
            return False

        parsed_nodes = []
        depth = 0
        for idx, m in enumerate(matches):
            if m.group(3):
                depth = max(0, depth - 1)
                continue

            raw_name = html.unescape(m.group(2)).strip()
            next_pos = matches[idx + 1].start() if idx + 1 < len(matches) else len(krona_xml)
            node_header = krona_xml[m.end():next_pos]

            val_m = val_re.search(node_header)
            hits_val = 0
            if val_m:
                try:
                    hits_val = int(round(float(val_m.group(1))))
                except ValueError:
                    hits_val = 0

            rank_m = rank_re.search(node_header)
            rank_str = html.unescape(rank_m.group(1)).strip() if rank_m else "-"

            tax_m = taxon_re.search(node_header)
            tax_str = html.unescape(tax_m.group(1)).strip() if tax_m else "-"

            parsed_nodes.append((depth, raw_name, rank_str, tax_str, hits_val))
            depth += 1

        if not parsed_nodes:
            return False

        self.txt_reader.insert("end", "=" * 84 + "\n", "h2")
        self.txt_reader.insert("end", "  Krona Interactive Taxonomic Composition Report\n", "h1")
        self.txt_reader.insert("end", "=" * 84 + "\n\n", "h2")

        max_mag = max((v for _, _, _, _, v in parsed_nodes), default=1) or 1
        table_rows = [["Taxon / Clade Hierarchy", "Rank", "TaxID", "Reads / Hits", "Relative Share", "Visual Proportion"]]
        for node_depth, taxon_name, rank_str, tax_str, m_val in parsed_nodes[:160]:
            indent = "  " * min(node_depth, 8)
            prefix = "|- " if node_depth > 0 else ""
            disp_name = f"{indent}{prefix}{taxon_name}"
            pct = (m_val / max_mag) * 100.0 if max_mag > 0 else 0.0
            filled = max(0, min(20, int(round((pct / 100.0) * 20))))
            bar_str = "#" * filled + "." * (20 - filled)
            table_rows.append([disp_name, rank_str, tax_str, f"{m_val:,}", f"{pct:5.1f}%", f"[{bar_str}]"])

        self._insert_formatted_box_table(
            f"Taxonomic Hierarchy ({len(parsed_nodes)} total clades, Root total = {max_mag:,} hits)",
            table_rows,
        )
        return True

    def _render_html_in_reader(self, file_path: Path, content: str):
        """
        Renders `.html` reports visually inside `self.txt_reader` with spacious word-wrapped
        sections, boxed tables, decoded MultiQC/Krona bar charts, and embedded PNG plots.
        """
        self.txt_reader.insert(
            "end",
            f"  [IN-APP HTML VISUAL VIEWER]  {file_path.name}   "
            f"(Click [Open in Browser / External App] above for full browser view)\n\n",
            "muted",
        )

        # 1. Check for Krona XML report
        if "krona" in file_path.name.lower() and self._render_krona_hierarchy(content):
            return

        # 2. Check for MultiQC compressed JSON report
        is_mqc = "multiqc" in file_path.name.lower() and self._render_multiqc_compressed_data(content)

        # 3. Extract <body>...</body> using fast O(1) string indexing instead of regex backtracking
        body_start_tag = content.find("<body")
        if body_start_tag == -1:
            body_start_tag = content.find("<BODY")
        if body_start_tag != -1:
            body_gt = content.find(">", body_start_tag)
            body_end = content.rfind("</body>")
            if body_end == -1:
                body_end = content.rfind("</BODY>")
            if body_gt != -1 and body_end > body_gt:
                body_html = content[body_gt + 1 : body_end]
            elif body_gt != -1:
                body_html = content[body_gt + 1 :]
            else:
                body_html = content
        else:
            body_html = content

        outer_self = self

        class _VisualHTMLDocumentParser(HTMLParser):
            def __init__(self):
                super().__init__()
                self.skip_depth = 0
                self.in_table = False
                self.in_caption = False
                self.caption_parts = []
                self.table_rows = []
                self.current_row = []
                self.in_cell = False
                self.current_cell = []
                self.heading_tag = None
                self.in_bold = False
                self.in_code = False
                self.in_list_item = False
                self.text_buf = []

            def _flush_text(self):
                if not self.text_buf:
                    return
                # Coalesce consecutive chunks with identical tag_name to reduce Tk insert() calls
                cur_tag = self.text_buf[0][1]
                cur_pieces = []
                for txt_piece, tag_name in self.text_buf:
                    if tag_name != cur_tag and cur_pieces:
                        if cur_tag:
                            outer_self.txt_reader.insert("end", "".join(cur_pieces), cur_tag)
                        else:
                            outer_self.txt_reader.insert("end", "".join(cur_pieces))
                        cur_pieces.clear()
                    cur_tag = tag_name
                    cur_pieces.append(txt_piece)
                if cur_pieces:
                    if cur_tag:
                        outer_self.txt_reader.insert("end", "".join(cur_pieces), cur_tag)
                    else:
                        outer_self.txt_reader.insert("end", "".join(cur_pieces))
                self.text_buf.clear()

            def handle_starttag(self, tag, attrs):
                t = tag.lower()
                attr_dict = dict(attrs)
                if t in ("script", "style", "svg", "noscript", "nav", "footer"):
                    self.skip_depth += 1
                    return
                if self.skip_depth > 0:
                    return

                # Skip empty floating TOC container in RMarkdown reports
                if t == "div" and attr_dict.get("id") == "TOC":
                    self.skip_depth += 1
                    return

                if t == "img":
                    src = attr_dict.get("src", "")
                    if "data:image/png;base64," in src:
                        self._flush_text()
                        outer_self._insert_base64_png_image(src)
                    return

                if t == "table":
                    self._flush_text()
                    self.in_table = True
                    self.table_rows = []
                    self.caption_parts = []
                elif self.in_table:
                    if t == "caption":
                        self.in_caption = True
                    elif t == "tr":
                        self.current_row = []
                    elif t in ("th", "td"):
                        self.in_cell = True
                        self.current_cell = []
                else:
                    if t == "h1":
                        self._flush_text()
                        self.heading_tag = "h1"
                        outer_self.txt_reader.insert("end", "\n" + "=" * 84 + "\n  ", "h2")
                    elif t == "h2":
                        self._flush_text()
                        self.heading_tag = "h1"
                        outer_self.txt_reader.insert("end", "\n" + "-" * 84 + "\n  ", "h2")
                    elif t in ("h3", "h4"):
                        self._flush_text()
                        self.heading_tag = "h3"
                        outer_self.txt_reader.insert("end", "\n  >> ", "h2")
                    elif t in ("strong", "b"):
                        self.in_bold = True
                    elif t in ("code", "pre"):
                        self.in_code = True
                    elif t == "li":
                        self._flush_text()
                        self.in_list_item = True
                        outer_self.txt_reader.insert("end", "\n    *  ", "key")
                    elif t in ("p", "br"):
                        self.text_buf.append(("\n  ", None))

            def handle_endtag(self, tag):
                t = tag.lower()
                if t in ("script", "style", "svg", "noscript", "nav", "footer") and self.skip_depth > 0:
                    self.skip_depth -= 1
                    return
                if self.skip_depth > 0:
                    if t == "div":
                        self.skip_depth -= 1
                    return

                if self.in_table:
                    if t == "caption":
                        self.in_caption = False
                    elif t in ("th", "td"):
                        self.in_cell = False
                        self.current_row.append(" ".join(self.current_cell).strip())
                    elif t == "tr":
                        if any(self.current_row):
                            self.table_rows.append(self.current_row)
                        self.current_row = []
                    elif t == "table":
                        self.in_table = False
                        cap_str = " ".join(self.caption_parts).strip()
                        outer_self._insert_formatted_box_table(cap_str, self.table_rows)
                        self.table_rows = []
                        self.caption_parts = []
                else:
                    if t == "h1":
                        self._flush_text()
                        outer_self.txt_reader.insert("end", "\n" + "=" * 84 + "\n", "h2")
                        self.heading_tag = None
                    elif t == "h2":
                        self._flush_text()
                        outer_self.txt_reader.insert("end", "\n" + "-" * 84 + "\n", "h2")
                        self.heading_tag = None
                    elif t in ("h3", "h4"):
                        self._flush_text()
                        outer_self.txt_reader.insert("end", "\n", None)
                        self.heading_tag = None
                    elif t in ("strong", "b"):
                        self.in_bold = False
                    elif t in ("code", "pre"):
                        self.in_code = False
                    elif t == "li":
                        self._flush_text()
                        self.in_list_item = False
                    elif t == "p":
                        self.text_buf.append(("\n", None))
                        self._flush_text()

            def handle_data(self, data):
                if self.skip_depth > 0:
                    return
                cleaned = " ".join(data.split())
                if not cleaned:
                    return
                if self.in_table:
                    if self.in_caption:
                        self.caption_parts.append(cleaned)
                    elif self.in_cell:
                        self.current_cell.append(cleaned)
                else:
                    if self.heading_tag:
                        self.text_buf.append((cleaned + " ", self.heading_tag))
                    elif self.in_code:
                        self.text_buf.append((cleaned + " ", "code"))
                    elif self.in_bold:
                        self.text_buf.append((cleaned + " ", "bold"))
                    else:
                        self.text_buf.append((cleaned + " ", None))

        if not is_mqc:
            parser = _VisualHTMLDocumentParser()
            try:
                parser.feed(body_html)
                parser._flush_text()
            except Exception:
                pass

    # =========================================================================
    # DARK MATTER (RdRp & RVDB) HIT-TO-FASTA LINKER & THRESHOLD EXTRACTOR
    # =========================================================================
    @staticmethod
    def _extract_contig_id_from_orf_label(label: str) -> str:
        """
        Extracts the original Contig ID from an EMBOSS getorf / CD-HIT ORF label
        (e.g. `gc_1_k79_10967_1_[2_-_2197]_flag=0...` -> `k79_10967`).
        """
        if not label:
            return "-"
        m = re.search(r"^gc_\d+_(.+?)_\d+_\[", label.strip())
        if m:
            return m.group(1)
        return label.strip().split()[0]

    @staticmethod
    def _parse_fasta_dict(
        fasta_paths: list[Path], target_ids: set[str] | None = None
    ) -> dict[str, tuple[str, str]]:
        """
        Reads one or more FASTA files and returns `{seq_id: (full_header, sequence_str)}`.
        If `target_ids` is provided, only stores matching sequence IDs and stops reading
        immediately once all `target_ids` have been found (avoiding scanning huge assembly FASTAs).
        """
        records = {}
        if target_ids is not None and len(target_ids) == 0:
            return records

        for fp in fasta_paths:
            if target_ids is not None and len(records) >= len(target_ids):
                break
            if not fp or not fp.is_file():
                continue
            try:
                with open(fp, "r", encoding="utf-8", errors="replace") as f:
                    cur_id = None
                    cur_hdr = ""
                    cur_seq = []
                    keep_cur = False
                    for raw_line in f:
                        if raw_line.startswith(">"):
                            if keep_cur and cur_id and cur_id not in records:
                                records[cur_id] = (cur_hdr, "".join(cur_seq))
                                if target_ids is not None and len(records) >= len(target_ids):
                                    keep_cur = False
                                    break
                            cur_hdr = raw_line[1:].strip()
                            cur_id = cur_hdr.split()[0] if cur_hdr else None
                            keep_cur = bool(cur_id) and (target_ids is None or cur_id in target_ids) and (cur_id not in records)
                            cur_seq = []
                        elif keep_cur:
                            cur_seq.append(raw_line.strip())
                    if keep_cur and cur_id and cur_id not in records:
                        records[cur_id] = (cur_hdr, "".join(cur_seq))
            except Exception:
                pass
        return records

    @staticmethod
    def _format_fasta_entry(header: str, seq: str, width: int = 80) -> str:
        clean_seq = re.sub(r"\s+", "", seq or "")
        wrapped = "\n".join(clean_seq[i : i + width] for i in range(0, len(clean_seq), width))
        return f">{header}\n{wrapped}\n"

    def _discover_darkmatter_contexts(self) -> list[dict]:
        """
        Finds all `(sample, tool)` Dark Matter / RVDB directories in `results/`.
        If the currently selected file in Tab 4 belongs to a specific sample, places that sample first.
        """
        out_dir = Path(self.var_output.get().strip() or str(PROJECT_DIR / "results"))
        if not out_dir.is_dir():
            return []

        contexts_map = {}
        for sample_dir in sorted(out_dir.iterdir()):
            if not sample_dir.is_dir() or sample_dir.name in ("benchmarks", "multiqc", "kraken2_all"):
                continue
            sample = sample_dir.name
            for sub in sorted(sample_dir.iterdir()):
                if not sub.is_dir():
                    continue
                if sub.name.startswith("darkmatter_"):
                    tool = sub.name[len("darkmatter_") :]
                elif sub.name.startswith("rvdb_"):
                    tool = sub.name[len("rvdb_") :]
                else:
                    continue
                key = (sample, tool)
                if key not in contexts_map:
                    dm_dir = sample_dir / f"darkmatter_{tool}"
                    rvdb_dir = sample_dir / f"rvdb_{tool}"
                    contexts_map[key] = {
                        "sample": sample,
                        "tool": tool,
                        "sample_dir": sample_dir,
                        "dm_dir": dm_dir,
                        "rvdb_dir": rvdb_dir,
                    }

        contexts = list(contexts_map.values())
        # If the user selected a file inside a specific sample/tool folder, prioritize that sample/tool
        if self.selected_result_path:
            sel_str = str(self.selected_result_path.resolve())
            prioritized = [
                c
                for c in contexts
                if str(c["dm_dir"].resolve()) in sel_str or str(c["rvdb_dir"].resolve()) in sel_str
            ]
            if not prioritized:
                prioritized = [c for c in contexts if str(c["sample_dir"].resolve()) in sel_str]
            if prioritized:
                return prioritized

        return contexts

    def _collect_linked_hits_for_context(
        self, ctx: dict, minscore: float, max_evalue: float, rvdb_cat: str
    ) -> dict:
        """
        Parses `{sample}_RdRp.tsv` and `{sample}_RVDB_Summary.csv` FIRST (small KB files),
        filters hits by `rdrp >= minscore` and `E_value <= max_evalue`, and then performs
        a targeted early-exit lookup for only the matching ORF and Contig FASTA sequences.
        """
        sample = ctx["sample"]
        tool = ctx["tool"]
        sample_dir = ctx["sample_dir"]
        dm_dir = ctx["dm_dir"]
        rvdb_dir = ctx["rvdb_dir"]

        target_orf_ids = set()
        target_contig_ids = set()
        inline_orf_seqs = {}

        # 1. Parse RdRp hits (`{sample}_fev2tsv_single_{tool}.tsv` or legacy `{sample}_RdRp.tsv`) first
        rdrp_tsv = next(
            (
                p
                for p in (
                    dm_dir / f"{sample}_fev2tsv_single_{tool}.tsv",
                    dm_dir / f"{sample}_{tool}_RdRp.tsv",
                    dm_dir / f"{sample}_RdRp.tsv",
                )
                if p.is_file()
            ),
            dm_dir / f"{sample}_fev2tsv_single_{tool}.tsv",
        )
        rdrp_hits = []
        total_rdrp_raw = 0
        if rdrp_tsv.is_file():
            try:
                lines = [ln for ln in rdrp_tsv.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
                if len(lines) > 1:
                    hdr = [h.strip() for h in lines[0].split("\t")]
                    col_map = {name: i for i, name in enumerate(hdr)}
                    for ln in lines[1:]:
                        total_rdrp_raw += 1
                        parts = ln.split("\t")
                        label = parts[col_map["Label"]].strip() if "Label" in col_map and col_map["Label"] < len(parts) else ""
                        raw_sc = parts[col_map["rdrp"]].strip() if "rdrp" in col_map and col_map["rdrp"] < len(parts) else "0"
                        pssm_sc = parts[col_map["pssm_score"]].strip() if "pssm_score" in col_map and col_map["pssm_score"] < len(parts) else "-"
                        gdd = parts[col_map["gdd"]].strip() if "gdd" in col_map and col_map["gdd"] < len(parts) else "-"
                        group = parts[col_map["pssm_group"]].strip() if "pssm_group" in col_map and col_map["pssm_group"] < len(parts) else "-"
                        aaseq = parts[col_map["aaseq"]].strip() if "aaseq" in col_map and col_map["aaseq"] < len(parts) else ""
                        try:
                            score_f = float(raw_sc)
                        except ValueError:
                            score_f = 0.0

                        if score_f >= minscore and label:
                            contig_id = self._extract_contig_id_from_orf_label(label)
                            target_orf_ids.add(label)
                            if contig_id != "-":
                                target_contig_ids.add(contig_id)
                            if aaseq:
                                inline_orf_seqs[label] = (f"{label} RdRp_score={score_f} motif_GDD={gdd}", aaseq)
                            rdrp_hits.append(
                                {
                                    "orf_id": label,
                                    "contig_id": contig_id,
                                    "rdrp_score": score_f,
                                    "pssm_score": pssm_sc,
                                    "gdd": gdd,
                                    "group": group,
                                }
                            )
            except Exception:
                pass

        rdrp_hits.sort(key=lambda x: -x["rdrp_score"])

        # 2. Parse RVDB hits (`{sample}_rvdb_summarize_{tool}.tsv` or legacy `{sample}_RVDB_Summary.csv`) first
        rvdb_csv = next(
            (
                p
                for p in (
                    rvdb_dir / f"{sample}_rvdb_summarize_{tool}.tsv",
                    rvdb_dir / f"{sample}_{tool}_RVDB_Summary.csv",
                    rvdb_dir / f"{sample}_RVDB_Summary.csv",
                )
                if p.is_file()
            ),
            rvdb_dir / f"{sample}_rvdb_summarize_{tool}.tsv",
        )
        rvdb_hits = []
        total_rvdb_raw = 0
        if rvdb_csv.is_file():
            try:
                lines = [ln for ln in rvdb_csv.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
                if len(lines) > 1:
                    hdr = [h.strip() for h in lines[0].split("\t")]
                    col_map = {name: i for i, name in enumerate(hdr)}
                    for ln in lines[1:]:
                        total_rvdb_raw += 1
                        parts = ln.split("\t")
                        seq_id = parts[col_map["Sequence_ID"]].strip() if "Sequence_ID" in col_map and col_map["Sequence_ID"] < len(parts) else ""
                        annot = parts[col_map["Annotation"]].strip() if "Annotation" in col_map and col_map["Annotation"] < len(parts) else "-"
                        sci_name = parts[col_map["Scientific_Name"]].strip() if "Scientific_Name" in col_map and col_map["Scientific_Name"] < len(parts) else "-"
                        conf = parts[col_map["Confidence"]].strip() if "Confidence" in col_map and col_map["Confidence"] < len(parts) else "-"
                        eval_str = parts[col_map["E_value"]].strip() if "E_value" in col_map and col_map["E_value"] < len(parts) else "1.0"
                        bit_str = parts[col_map["Bit_Score"]].strip() if "Bit_Score" in col_map and col_map["Bit_Score"] < len(parts) else "0"
                        hmm_fam = parts[col_map["HMM_Family"]].strip() if "HMM_Family" in col_map and col_map["HMM_Family"] < len(parts) else "-"
                        try:
                            eval_f = float(eval_str)
                        except ValueError:
                            eval_f = 1.0

                        is_pol = any(
                            kw in annot.lower()
                            for kw in ("rdrp", "polymerase", "replicase", "reverse transcriptase", "rna-dependent")
                        )
                        is_struct = any(
                            kw in annot.lower()
                            for kw in ("capsid", "coat", "envelope", "glycoprotein", "nucleocapsid", "matrix", "structural", "virion")
                        )
                        cat_label = "Polymerase/RdRp" if is_pol else ("Structural" if is_struct else "Other Viral")

                        if rvdb_cat == "Polymerase / RdRp" and not is_pol:
                            continue
                        if rvdb_cat == "Structural Only" and not is_struct:
                            continue

                        if eval_f <= max_evalue and seq_id:
                            contig_id = self._extract_contig_id_from_orf_label(seq_id)
                            target_orf_ids.add(seq_id)
                            if contig_id != "-":
                                target_contig_ids.add(contig_id)
                            rvdb_hits.append(
                                {
                                    "orf_id": seq_id,
                                    "contig_id": contig_id,
                                    "category": cat_label,
                                    "confidence": conf,
                                    "evalue": eval_f,
                                    "evalue_str": eval_str,
                                    "bit_score": bit_str,
                                    "hmm_family": hmm_fam,
                                    "scientific_name": sci_name,
                                    "annotation": annot,
                                }
                            )
            except Exception:
                pass

        rvdb_hits.sort(key=lambda x: x["evalue"])

        # 3. Targeted lookup in ORF FASTA sources (new {sample}_{rule} + legacy names, small files FIRST!)
        orf_sources = [
            dm_dir / f"{sample}_dm_validate_{tool}_orfs.fasta",
            dm_dir / f"{sample}_palm_annot_{tool}.fasta",
            rvdb_dir / f"{sample}_rvdb_validate_pol_{tool}_orfs.fasta",
            rvdb_dir / f"{sample}_rvdb_validate_structural_{tool}_orfs.fasta",
            dm_dir / f"{sample}_cd_hit_{tool}.fasta",
            dm_dir / f"{sample}_RdRp_Orfs.fasta",
            dm_dir / f"{sample}_RdRp.fasta",
            rvdb_dir / f"{sample}_pol_Orfs.fasta",
            rvdb_dir / f"{sample}_structural_Orfs.fasta",
            dm_dir / f"{sample}_ORFs.fasta",
        ]
        orf_seqs = self._parse_fasta_dict(orf_sources, target_ids=target_orf_ids)
        for oid, val in inline_orf_seqs.items():
            if oid not in orf_seqs:
                orf_seqs[oid] = val

        # 4. Targeted lookup in Contig FASTA sources (new {sample}_{rule} + legacy names, small files FIRST!)
        contig_sources = [
            dm_dir / f"{sample}_dm_validate_{tool}_contigs.fasta",
            rvdb_dir / f"{sample}_rvdb_validate_pol_{tool}_contigs.fasta",
            rvdb_dir / f"{sample}_rvdb_validate_structural_{tool}_contigs.fasta",
            dm_dir / f"{sample}_get_nohit_fasta_{tool}.fasta",
            sample_dir / tool / f"{sample}_{tool}.fasta",
            dm_dir / f"{sample}_RdRp_contigs.fasta",
            rvdb_dir / f"{sample}_pol_contigs.fasta",
            rvdb_dir / f"{sample}_structural_contigs.fasta",
            dm_dir / f"{sample}_{tool}_nohit.fasta",
            sample_dir / tool / "final.contigs.fa",
            sample_dir / tool / "contigs.fasta",
            sample_dir / tool / "assembly.fasta",
            sample_dir / tool / "consensus.fasta",
        ]
        if tool.startswith("spades_k"):
            kval = tool[len("spades_k") :]
            contig_sources.insert(4, sample_dir / "spades" / f"kmer_{kval}" / f"{sample}_spades_k{kval}.fasta")
            contig_sources.append(sample_dir / "spades" / f"kmer_{kval}" / "contigs.fasta")
        elif tool.startswith("medaka_"):
            asm_base = tool.split("medaka_", 1)[1]
            contig_sources.insert(4, sample_dir / tool / f"{sample}_medaka_polish_{asm_base}.fasta")
        contig_seqs = self._parse_fasta_dict(contig_sources, target_ids=target_contig_ids)

        for h in rdrp_hits:
            h["has_orf_seq"] = h["orf_id"] in orf_seqs
            h["has_contig_seq"] = h["contig_id"] in contig_seqs
        for h in rvdb_hits:
            h["has_orf_seq"] = h["orf_id"] in orf_seqs
            h["has_contig_seq"] = h["contig_id"] in contig_seqs

        return {
            "sample": sample,
            "tool": tool,
            "dm_dir": dm_dir,
            "rvdb_dir": rvdb_dir,
            "total_rdrp_raw": total_rdrp_raw,
            "rdrp_hits": rdrp_hits,
            "total_rvdb_raw": total_rvdb_raw,
            "rvdb_hits": rvdb_hits,
            "orf_seqs": orf_seqs,
            "contig_seqs": contig_seqs,
        }

    def preview_linked_darkmatter_hits(self):
        """
        Links RdRp and RVDB hits above the user-selected threshold to their ORF and Contig
        FASTA sequences and renders an interactive summary + FASTA viewer inside `self.txt_reader`.
        """
        try:
            minscore = float(self.var_dm_filter_minscore.get())
        except (ValueError, tk.TclError):
            minscore = 75.0
        try:
            max_evalue = float(self.var_rvdb_max_evalue.get())
        except (ValueError, tk.TclError):
            max_evalue = 1e-5
        rvdb_cat = self.var_rvdb_cat.get()

        contexts = self._discover_darkmatter_contexts()
        self.txt_reader.configure(wrap="none", padx=14, pady=12, spacing1=2, spacing2=1, spacing3=2)
        self.txt_reader.delete("1.0", "end")
        self.var_reader_badge.set("[LINKED RdRp & RVDB FASTAs]")
        self.var_reader_title.set(
            f"Dark Matter Hit-to-FASTA Linker (RdRp >= {minscore:g} | RVDB E <= {max_evalue:g})"
        )

        if not contexts:
            self.txt_reader.insert(
                "end",
                "No Dark Matter or RVDB output directories found in results/ yet.\n"
                "Run the pipeline with --darkmatter (and optionally --rvdb) first.\n",
                "error",
            )
            return

        self.txt_reader.insert("end", "=" * 96 + "\n", "h2")
        self.txt_reader.insert(
            "end",
            f"  DARK MATTER (RdRp & RVDB) HIT-TO-FASTA LINKER & THRESHOLD VIEWER\n"
            f"  Active Filter: RdRp MinScore >= {minscore:g}  |  RVDB E-value <= {max_evalue:g} ({rvdb_cat})\n",
            "h1",
        )
        self.txt_reader.insert("end", "=" * 96 + "\n\n", "h2")

        for ctx in contexts:
            res = self._collect_linked_hits_for_context(ctx, minscore, max_evalue, rvdb_cat)
            sample = res["sample"]
            tool = res["tool"]
            rdrp_hits = res["rdrp_hits"]
            rvdb_hits = res["rvdb_hits"]
            orf_seqs = res["orf_seqs"]
            contig_seqs = res["contig_seqs"]

            rdrp_contigs = {h["contig_id"] for h in rdrp_hits if h["contig_id"] != "-"}
            rvdb_contigs = {h["contig_id"] for h in rvdb_hits if h["contig_id"] != "-"}
            dual_contigs = sorted(rdrp_contigs & rvdb_contigs)

            self.txt_reader.insert("end", f"  SAMPLE: {sample}   |   ASSEMBLER: {tool}\n", "h1")
            self.txt_reader.insert(
                "end",
                f"  RdRp Hits >= {minscore:g}: {len(rdrp_hits)}/{res['total_rdrp_raw']} ORFs ({len(rdrp_contigs)} Contigs)   |   "
                f"RVDB Hits <= {max_evalue:g}: {len(rvdb_hits)}/{res['total_rvdb_raw']} ORFs ({len(rvdb_contigs)} Contigs)\n\n",
                "bar_fill",
            )

            # 1. RdRp Hits Table
            if rdrp_hits:
                r_rows = [
                    [
                        "ORF Label",
                        "Linked Contig ID",
                        "RdRp Score",
                        "PSSM Score",
                        "GDD Motif",
                        "Taxon / Group",
                        "FASTA Ready",
                    ]
                ]
                for h in rdrp_hits:
                    fasta_status = (
                        "ORF + Contig"
                        if (h["has_orf_seq"] and h["has_contig_seq"])
                        else ("ORF Only" if h["has_orf_seq"] else "Contig Only")
                    )
                    r_rows.append(
                        [
                            h["orf_id"],
                            h["contig_id"],
                            f"{h['rdrp_score']:.1f}",
                            str(h["pssm_score"]),
                            str(h["gdd"]),
                            str(h["group"]),
                            fasta_status,
                        ]
                    )
                self._insert_formatted_box_table(
                    f"{sample} ({tool}) - Linked RdRp Hits (Score >= {minscore:g})", r_rows
                )
            else:
                self.txt_reader.insert(
                    "end",
                    f"  [No RdRp hits with score >= {minscore:g} in {sample}/{tool}. "
                    f"Total candidates in TSV: {res['total_rdrp_raw']}. Try lowering RdRp MinScore above!]\n\n",
                    "muted",
                )

            # 2. RVDB Hits Table
            if rvdb_hits:
                v_rows = [
                    [
                        "ORF Sequence_ID",
                        "Linked Contig ID",
                        "Category",
                        "Confidence",
                        "E-value",
                        "BitScore",
                        "Scientific Name",
                        "FASTA Ready",
                    ]
                ]
                for h in rvdb_hits[:150]:
                    fasta_status = (
                        "ORF + Contig"
                        if (h["has_orf_seq"] and h["has_contig_seq"])
                        else ("ORF Only" if h["has_orf_seq"] else "Contig Only")
                    )
                    v_rows.append(
                        [
                            h["orf_id"],
                            h["contig_id"],
                            h["category"],
                            h["confidence"],
                            h["evalue_str"],
                            h["bit_score"],
                            h["scientific_name"],
                            fasta_status,
                        ]
                    )
                self._insert_formatted_box_table(
                    f"{sample} ({tool}) - Linked RVDB Hits (E-value <= {max_evalue:g})", v_rows
                )
            elif res["total_rvdb_raw"] > 0:
                self.txt_reader.insert(
                    "end",
                    f"  [No RVDB hits matching E-value <= {max_evalue:g} ({rvdb_cat}) in {sample}/{tool}.]\n\n",
                    "muted",
                )

            # 3. Dual-Validated Contigs (Confirmed by BOTH PalmAnnot RdRp & RVDB HMMs)
            if dual_contigs:
                self.txt_reader.insert(
                    "end",
                    f"  >> HIGH-CONFIDENCE CONTIGS MATCHED BY BOTH RdRp AND RVDB ({len(dual_contigs)} contigs): "
                    + ", ".join(dual_contigs)
                    + "\n\n",
                    "table_hdr",
                )

            # 4. Preview Linked Contig FASTA & ORF FASTA Sequences
            linked_orf_ids = []
            for h in rdrp_hits + rvdb_hits:
                if h["orf_id"] not in linked_orf_ids and h["orf_id"] in orf_seqs:
                    linked_orf_ids.append(h["orf_id"])

            linked_contig_ids = []
            for h in rdrp_hits + rvdb_hits:
                if h["contig_id"] not in linked_contig_ids and h["contig_id"] in contig_seqs:
                    linked_contig_ids.append(h["contig_id"])

            if linked_contig_ids:
                self.txt_reader.insert(
                    "end",
                    f"  LINKED CONTIG NUCLEOTIDE FASTA SEQUENCES ({len(linked_contig_ids)} contigs >= threshold):\n",
                    "h2",
                )
                self.txt_reader.insert("end", "  " + "-" * 88 + "\n", "table_border")
                for cid in linked_contig_ids[:50]:
                    hdr, seq = contig_seqs[cid]
                    self.txt_reader.insert("end", f">{hdr}  [Linked_Contig_Len={len(seq)}bp]\n", "fasta_hdr")
                    wrapped_seq = "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n"
                    self.txt_reader.insert("end", wrapped_seq)
                self.txt_reader.insert("end", "\n")

            if linked_orf_ids:
                self.txt_reader.insert(
                    "end",
                    f"  LINKED ORF AMINO ACID FASTA SEQUENCES ({len(linked_orf_ids)} ORFs >= threshold):\n",
                    "h2",
                )
                self.txt_reader.insert("end", "  " + "-" * 88 + "\n", "table_border")
                for oid in linked_orf_ids[:50]:
                    hdr, seq = orf_seqs[oid]
                    self.txt_reader.insert("end", f">{hdr}  [Linked_ORF_Len={len(seq)}aa]\n", "fasta_hdr")
                    wrapped_seq = "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n"
                    self.txt_reader.insert("end", wrapped_seq)
                self.txt_reader.insert("end", "\n")

    def extract_and_save_linked_fastas(self):
        """
        Extracts all linked ORF and Contig FASTA sequences for RdRp hits (`rdrp >= minscore`)
        and RVDB hits (`E_value <= max_evalue`), saves them into `results/<sample>/darkmatter_<tool>/`
        (and `rvdb_<tool>/`), refreshes the Results file list, and displays the extracted sequences.
        """
        try:
            minscore = float(self.var_dm_filter_minscore.get())
        except (ValueError, tk.TclError):
            minscore = 75.0
        try:
            max_evalue = float(self.var_rvdb_max_evalue.get())
        except (ValueError, tk.TclError):
            max_evalue = 1e-5
        rvdb_cat = self.var_rvdb_cat.get()

        contexts = self._discover_darkmatter_contexts()
        if not contexts:
            messagebox.showwarning(
                "No Dark Matter Outputs Found",
                "Could not find any `darkmatter_<tool>` or `rvdb_<tool>` folders in the Output Directory.",
            )
            return

        saved_files = []
        total_orfs_saved = 0
        total_contigs_saved = 0

        for ctx in contexts:
            res = self._collect_linked_hits_for_context(ctx, minscore, max_evalue, rvdb_cat)
            sample = res["sample"]
            dm_dir = res["dm_dir"]
            rvdb_dir = res["rvdb_dir"]
            rdrp_hits = res["rdrp_hits"]
            rvdb_hits = res["rvdb_hits"]
            orf_seqs = res["orf_seqs"]
            contig_seqs = res["contig_seqs"]

            # 1. Save Filtered RdRp ORFs & Contigs
            if dm_dir.is_dir():
                ms_tag = f"{int(minscore)}" if minscore.is_integer() else f"{minscore}"
                out_rdrp_orfs = dm_dir / f"{sample}_RdRp_minscore{ms_tag}_Orfs.fasta"
                out_rdrp_contigs = dm_dir / f"{sample}_RdRp_minscore{ms_tag}_contigs.fasta"
                std_rdrp_orfs = dm_dir / f"{sample}_RdRp_Orfs.fasta"
                std_rdrp_contigs = dm_dir / f"{sample}_RdRp_contigs.fasta"

                written_orfs = set()
                orf_fasta_chunks = []
                for h in rdrp_hits:
                    oid = h["orf_id"]
                    if oid not in written_orfs and oid in orf_seqs:
                        written_orfs.add(oid)
                        hdr, seq = orf_seqs[oid]
                        orf_fasta_chunks.append(
                            self._format_fasta_entry(
                                f"{hdr} | RdRp_score={h['rdrp_score']:.1f} | Contig={h['contig_id']}", seq
                            )
                        )

                written_contigs = set()
                contig_fasta_chunks = []
                for h in rdrp_hits:
                    cid = h["contig_id"]
                    if cid not in written_contigs and cid in contig_seqs:
                        written_contigs.add(cid)
                        hdr, seq = contig_seqs[cid]
                        contig_fasta_chunks.append(
                            self._format_fasta_entry(
                                f"{hdr} | Linked_RdRp_ORF={h['orf_id']} | RdRp_score={h['rdrp_score']:.1f}", seq
                            )
                        )

                orf_text = "".join(orf_fasta_chunks)
                contig_text = "".join(contig_fasta_chunks)

                out_rdrp_orfs.write_text(orf_text, encoding="utf-8", newline="\n")
                out_rdrp_contigs.write_text(contig_text, encoding="utf-8", newline="\n")
                std_rdrp_orfs.write_text(orf_text, encoding="utf-8", newline="\n")
                std_rdrp_contigs.write_text(contig_text, encoding="utf-8", newline="\n")

                saved_files.extend([out_rdrp_orfs, out_rdrp_contigs])
                total_orfs_saved += len(written_orfs)
                total_contigs_saved += len(written_contigs)

            # 2. Save Filtered RVDB ORFs & Contigs (if RVDB directory exists and has hits)
            if rvdb_dir.is_dir() and rvdb_hits:
                out_rvdb_orfs = rvdb_dir / f"{sample}_RVDB_filtered_Orfs.fasta"
                out_rvdb_contigs = rvdb_dir / f"{sample}_RVDB_filtered_contigs.fasta"

                v_orfs = set()
                v_orf_chunks = []
                for h in rvdb_hits:
                    oid = h["orf_id"]
                    if oid not in v_orfs and oid in orf_seqs:
                        v_orfs.add(oid)
                        hdr, seq = orf_seqs[oid]
                        v_orf_chunks.append(
                            self._format_fasta_entry(
                                f"{hdr} | RVDB_eval={h['evalue_str']} | {h['category']} | {h['scientific_name']}",
                                seq,
                            )
                        )

                v_contigs = set()
                v_contig_chunks = []
                for h in rvdb_hits:
                    cid = h["contig_id"]
                    if cid not in v_contigs and cid in contig_seqs:
                        v_contigs.add(cid)
                        hdr, seq = contig_seqs[cid]
                        v_contig_chunks.append(
                            self._format_fasta_entry(
                                f"{hdr} | Linked_RVDB_ORF={h['orf_id']} | RVDB_eval={h['evalue_str']} | {h['scientific_name']}",
                                seq,
                            )
                        )

                out_rvdb_orfs.write_text("".join(v_orf_chunks), encoding="utf-8", newline="\n")
                out_rvdb_contigs.write_text("".join(v_contig_chunks), encoding="utf-8", newline="\n")
                saved_files.extend([out_rvdb_orfs, out_rvdb_contigs])
                total_orfs_saved += len(v_orfs)
                total_contigs_saved += len(v_contigs)

        # Refresh file list and preview the linked hits & sequences
        self.refresh_results_list()
        self.preview_linked_darkmatter_hits()

        files_summary = "\n".join(f"- {p.name}" for p in saved_files[:8])
        messagebox.showinfo(
            "Linked FASTAs Extracted",
            f"Successfully extracted {total_orfs_saved} ORF sequence(s) and {total_contigs_saved} Contig sequence(s) "
            f"above threshold (RdRp >= {minscore:g}, RVDB E <= {max_evalue:g}).\n\n"
            f"Saved files:\n{files_summary}",
        )

    # =========================================================================
    # EXECUTION LOGIC & PROCESS MANAGEMENT
    # =========================================================================
    def on_mode_change(self):
        mode = self.mode_var.get()
        if mode == "slurm":
            self.var_profile.set("profile_slurm")
            self.load_profile_into_editor("profile_slurm", silent=True)
            self.append_terminal("[MODE CHANGED] Configured for Slurm Cluster (pipeline/profiles/profile_slurm).\n")
        else:
            self.var_profile.set("local")
            self.load_profile_into_editor("local", silent=True)
            self.append_terminal("[MODE CHANGED] Configured for Local Execution (pipeline/profiles/local).\n")

    def build_command_args(self, dry_run=False):
        args = ["DeepVir.py"]

        # Paths
        args += ["--input", self.var_input.get().strip()]
        args += ["--output", self.var_output.get().strip()]

        if self.var_temp.get().strip():
            args += ["--temp-dir", self.var_temp.get().strip()]

        if self.var_profile.get().strip():
            args += ["--profile", self.var_profile.get().strip()]

        args += ["--jobs", str(self.var_jobs.get())]

        # SRA File (if specified)
        if self.var_sra.get().strip():
            args += ["--sra", self.var_sra.get().strip()]

        # Direct SRA Accessions typed in the GUI text box
        direct_sra_tokens = self.get_direct_sra_tokens()
        if direct_sra_tokens:
            LOGS_DIR.mkdir(parents=True, exist_ok=True)
            sra_gui_file = LOGS_DIR / "sra_accessions_gui.txt"
            sra_gui_file.write_text("\n".join(direct_sra_tokens) + "\n", encoding="utf-8", newline="\n")
            args += ["--sra-list", ",".join(direct_sra_tokens)]

        # External Databases
        if self.var_diamond_db.get().strip():
            args += ["--diamond_db", self.var_diamond_db.get().strip()]
        if self.var_kraken2_db.get().strip():
            args += ["--kraken2_db", self.var_kraken2_db.get().strip()]

        # Read Preprocessing (fastp QC & Trimming)
        try:
            fastp_len_val = int(self.var_fastp_len.get())
            if fastp_len_val > 0:
                args += ["--fastp-length", str(fastp_len_val)]
        except (ValueError, tk.TclError):
            pass

        try:
            fastp_qual_val = int(self.var_fastp_qual.get())
            if fastp_qual_val > 0:
                args += ["--fastp-quality", str(fastp_qual_val)]
        except (ValueError, tk.TclError):
            pass

        # Pipeline Analysis Modules
        if self.mod_host_filter.get():
            args.append("--host-filter")
            host_tokens = self.get_host_accession_tokens()
            host_acc = ",".join(host_tokens) if host_tokens else "GCA_022682495.2"
            args += ["--host-accession", host_acc]

        if self.mod_assembly.get():
            args.append("--assembly")
            selected_assemblers = []
            if self.asm_spades.get():
                selected_assemblers.append("spades")
            if self.asm_megahit.get():
                selected_assemblers.append("megahit")
            if self.asm_flye.get():
                selected_assemblers.append("flye")
            if self.asm_raven.get():
                selected_assemblers.append("raven")
            if self.asm_medaka_flye.get():
                selected_assemblers.append("medaka_flye")
            if self.asm_medaka_raven.get():
                selected_assemblers.append("medaka_raven")

            if selected_assemblers:
                args += ["--assembler", ",".join(selected_assemblers)]

        if self.mod_reads_kraken.get():
            args.append("--reads-kraken")
        if self.mod_reads_diamond.get():
            args.append("--reads-diamond")
        if self.mod_kraken2.get():
            args.append("--kraken2")
        if self.mod_diamond.get():
            args.append("--diamond")
        if self.mod_darkmatter.get():
            args.append("--darkmatter")
            try:
                ms_val = int(self.var_rdrp_minscore.get())
                args += ["--rdrp-minscore", str(max(0, min(100, ms_val)))]
            except (ValueError, tk.TclError):
                args += ["--rdrp-minscore", "75"]
        if self.mod_rvdb.get():
            args.append("--rvdb")
        if self.mod_download_only.get():
            args.append("--download-only")
        if self.mod_remove_download.get():
            args.append("--remove-download")
        if dry_run:
            args.append("--dry-run")

        return args

    def start_execution(self, dry_run=False):
        if self.is_running:
            messagebox.showwarning("Warning", "An execution is already currently running!")
            return

        if self.mod_host_filter.get() and not self.var_host_accession.get().strip():
            messagebox.showwarning(
                "Host Accession Missing",
                "You enabled Host Genome Filtering (--host-filter) but did not specify an NCBI Host Accession.\n\n"
                "Please enter a valid NCBI genome assembly accession (e.g., GCA_022682495.2 for Desmodus rotundus).",
            )
            return

        if self.mod_assembly.get():
            selected_assemblers = [
                name
                for name, var in [
                    ("spades", self.asm_spades),
                    ("megahit", self.asm_megahit),
                    ("flye", self.asm_flye),
                    ("raven", self.asm_raven),
                    ("medaka_flye", self.asm_medaka_flye),
                    ("medaka_raven", self.asm_medaka_raven),
                ]
                if var.get()
            ]
            if not selected_assemblers:
                messagebox.showwarning(
                    "Assembler Missing",
                    "You enabled De Novo Assembly (--assembly) but did not select any assembler tools.\n\n"
                    "Please select at least one assembler engine (e.g., SPAdes, MEGAHIT, or Flye) "
                    "before starting the execution.",
                )
                return

        # Automatically persist any per-rule resource changes to the active profile before running
        self.save_editor_to_profile(show_popup=False)

        out_dir = Path(self.var_output.get().strip())
        out_dir.mkdir(parents=True, exist_ok=True)
        LOGS_DIR.mkdir(parents=True, exist_ok=True)

        log_file = LOGS_DIR / "deepvir_desktop_gui.log"
        self.active_log_path = log_file

        args = self.build_command_args(dry_run)

        self.notebook.select(self.tab_terminal)

        mode_str = "DRY-RUN (Preview)" if dry_run else "FULL WORKFLOW EXECUTION"
        self.append_terminal("=" * 70 + "\n")
        self.append_terminal(f" DeepVir Desktop GUI - {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.append_terminal(f" Execution Mode: {mode_str}\n")
        self.append_terminal(f" Command: python3 {' '.join(args)}\n")
        self.append_terminal(f" Log File: {log_file}\n")
        self.append_terminal("=" * 70 + "\n\n")

        self.set_running_state(True, dry_run)

        worker = threading.Thread(target=self.run_process_thread, args=(args, log_file, dry_run), daemon=True)
        worker.start()

    def run_process_thread(self, args, log_file, dry_run):
        try:
            executor_script = str(PROJECT_DIR / "DeepVir.py")
            cmd = [sys.executable, "-u", executor_script] + args[1:]

            with open(log_file, "a", encoding="utf-8", buffering=1) as lf:
                self.current_process = subprocess.Popen(
                    cmd,
                    cwd=str(PROJECT_DIR),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )

                for line in iter(self.current_process.stdout.readline, ""):
                    clean_disk_line = ANSI_ESCAPE_RE.sub("", line)
                    lf.write(clean_disk_line)
                    self.log_queue.put(line)

                self.current_process.stdout.close()
                return_code = self.current_process.wait()

            if return_code == 0:
                self.log_queue.put("\n" + "=" * 70 + "\n [OK] Pipeline finished with SUCCESS!\n" + "=" * 70 + "\n")
                self.after(0, lambda: self.on_execution_finished(True))
            else:
                self.log_queue.put(
                    f"\n" + "=" * 70 + f"\n [ERROR] Execution failed (Exit code: {return_code})\n" + "=" * 70 + "\n"
                )
                self.after(0, lambda: self.on_execution_finished(False))

        except Exception as e:
            self.log_queue.put(f"\n[EXECUTION ERROR]: {str(e)}\n")
            self.after(0, lambda: self.on_execution_finished(False))

    def on_execution_finished(self, success):
        self.set_running_state(False)
        self.refresh_results_list()
        if success:
            self.status_var.set("[ FINISHED ]")
            self.stage_var.set("Stage: Completed")
            self.lbl_status_badge.config(fg="#10b981", bg="#064e3b")
            messagebox.showinfo("Success", "DeepVir pipeline finished successfully!")
        else:
            self.status_var.set("[ FAILED ]")
            self.stage_var.set("Stage: Failed")
            self.lbl_status_badge.config(fg="#ef4444", bg="#450a0a")

    def stop_execution(self):
        if self.current_process and self.is_running:
            try:
                self.current_process.terminate()
                self.append_terminal("\n[WARNING] Termination signal sent to process.\n")
                self.status_var.set("[ ABORTED ]")
                self.stage_var.set("Stage: Aborted")
                self.lbl_status_badge.config(fg="#ef4444", bg="#450a0a")
            except Exception as e:
                self.append_terminal(f"\n[ERROR STOPPING PROCESS]: {e}\n")
            finally:
                self.set_running_state(False)

    def set_running_state(self, is_running, dry_run=False):
        self.is_running = is_running
        if is_running:
            text = "[ RUNNING DRY-RUN ]" if dry_run else "[ RUNNING PIPELINE ]"
            self.status_var.set(text)
            self.stage_var.set("Stage: Initializing...")
            self.lbl_status_badge.config(fg="#f59e0b", bg="#451a03")
            self.btn_run.config(state="disabled")
            self.btn_dryrun.config(state="disabled")
            self.btn_stop.config(state="normal")
        else:
            self.btn_run.config(state="normal")
            self.btn_dryrun.config(state="normal")
            self.btn_stop.config(state="disabled")

    def process_log_queue(self):
        while not self.log_queue.empty():
            line = self.log_queue.get_nowait()
            self.append_terminal(line)
        self.after(100, self.process_log_queue)

    def update_stage_from_line(self, clean_line):
        stripped = clean_line.strip()
        if not stripped:
            return
        if "Checking Conda environment" in stripped:
            self.stage_var.set("Stage: Checking Conda Env")
        elif "Generating sample list" in stripped:
            self.stage_var.set("Stage: Scanning Samples")
        elif "Unlocking working directory" in stripped:
            self.stage_var.set("Stage: Unlocking Workdir")
        elif "Performing a dry-run" in stripped:
            self.stage_var.set("Stage: Snakemake Dry-Run (DAG)")
        elif "Starting main execution" in stripped:
            self.stage_var.set("Stage: Main Execution")
        elif stripped.startswith("Creating conda environment"):
            env_name = stripped.split("/")[-1].replace(".yaml...", "")
            self.stage_var.set(f"Stage: Creating Conda Env ({env_name})")
        elif stripped.startswith("rule ") and stripped.endswith(":"):
            rule_name = stripped[5:-1]
            self.stage_var.set(f"Rule: {rule_name}")
        elif "steps (" in stripped and ") done" in stripped:
            self.stage_var.set(f"Progress: {stripped}")
        elif "Creating DAG & Report" in stripped:
            self.stage_var.set("Stage: Generating HTML Report")

    def append_terminal(self, text):
        clean_text = ANSI_ESCAPE_RE.sub("", text)
        self.update_stage_from_line(clean_text)

        tag = None
        if "[ERROR]" in clean_text or "Execution failed" in clean_text or "FAILED" in clean_text:
            tag = "error"
        elif "[WARNING]" in clean_text:
            tag = "warn"
        elif "[INFO]" in clean_text:
            tag = "info"
        elif clean_text.strip().startswith("rule ") and clean_text.strip().endswith(":"):
            tag = "rule"
        elif (
            clean_text.strip().startswith("> Snakemake:")
            or "[OK]" in clean_text
            or "OK (" in clean_text
            or "Active" in clean_text
        ):
            tag = "success"

        if tag:
            self.txt_terminal.insert("end", clean_text, tag)
        else:
            self.txt_terminal.insert("end", clean_text)

        if self.var_autoscroll.get():
            self.txt_terminal.see("end")

        lines = int(self.txt_terminal.index("end-1c").split(".")[0])
        self.lbl_lines.config(text=f"{lines} lines")

    def clear_terminal(self):
        self.txt_terminal.delete("1.0", "end")
        self.lbl_lines.config(text="0 lines")
        self.stage_var.set("Stage: Idle")

    def save_log_file(self):
        LOGS_DIR.mkdir(parents=True, exist_ok=True)
        dest = filedialog.asksaveasfilename(
            initialdir=str(LOGS_DIR),
            defaultextension=".log",
            filetypes=[("Log Files", "*.log"), ("Text Files", "*.txt")],
        )
        if dest:
            with open(dest, "w", encoding="utf-8") as f:
                f.write(self.txt_terminal.get("1.0", "end-1c"))
            messagebox.showinfo("Saved", f"Log saved successfully to:\n{dest}")

    # =========================================================================
    # OUTPUT LAUNCHERS & SYSTEM INTEGRATION (WSL, WINDOWS, MACOS, LINUX)
    # =========================================================================
    def _is_wsl(self) -> bool:
        if os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP"):
            return True
        try:
            return "microsoft" in platform.uname().release.lower()
        except Exception:
            return False

    def _to_windows_path(self, target_path: Path) -> str:
        """Converts a WSL POSIX path (e.g. `/mnt/d/DeepVir/results/...`) to a Windows path (`D:\\DeepVir\\results\\...`)."""
        posix_str = str(target_path.resolve()).replace("\\", "/")
        m = re.match(r"^/mnt/([a-zA-Z])/(.*)$", posix_str)
        if m:
            drive = m.group(1).upper()
            rest = m.group(2).replace("/", "\\")
            return f"{drive}:\\{rest}"
        try:
            out = subprocess.check_output(["wslpath", "-w", str(target_path.resolve())], text=True).strip()
            if out:
                return out
        except Exception:
            pass
        return str(target_path)

    def open_path_in_os(self, target_path: Path):
        """Opens any file (HTML report, TSV, TXT, etc.) in the host OS default browser or viewer."""
        if not target_path or not target_path.exists():
            messagebox.showwarning("Not Found", f"File not found:\n{target_path}")
            return

        sys_name = platform.system()
        if sys_name == "Windows":
            os.startfile(str(target_path))
        elif self._is_wsl():
            win_path = self._to_windows_path(target_path)
            subprocess.Popen(["cmd.exe", "/c", "start", "", win_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elif sys_name == "Darwin":
            subprocess.Popen(["open", str(target_path)])
        else:
            try:
                subprocess.Popen(["xdg-open", str(target_path)])
            except Exception:
                webbrowser.open(target_path.resolve().as_uri())

    def open_selected_in_external_app(self):
        """Opens the currently selected file in Tab 4 in the host OS browser or default application."""
        if not self.selected_result_path or not self.selected_result_path.exists():
            messagebox.showwarning("No File Selected", "Please select an output file from the list on the left first.")
            return
        self.open_path_in_os(self.selected_result_path)

    def open_output_file(self, filename):
        out_dir = Path(self.var_output.get().strip())
        matches = list(out_dir.glob(f"**/{filename}"))
        if matches:
            self.load_file_into_reader(matches[0])
            self.open_path_in_os(matches[0])
        else:
            messagebox.showwarning(
                "Not Found",
                f"The file '{filename}' was not found in the output directory yet.\n"
                "Please verify that the analysis has completed.",
            )

    def open_output_pattern(self, pattern):
        out_dir = Path(self.var_output.get().strip())
        matches = list(out_dir.glob(f"**/{pattern}"))
        if matches:
            self.load_file_into_reader(matches[0])
            self.open_path_in_os(matches[0])
        else:
            messagebox.showwarning("Not Found", f"No files matching pattern '{pattern}' were found.")

    def open_folder_in_os(self, folder_path: Path):
        folder_path.mkdir(parents=True, exist_ok=True)
        sys_name = platform.system()
        if sys_name == "Windows":
            os.startfile(str(folder_path))
        elif self._is_wsl():
            win_path = self._to_windows_path(folder_path)
            subprocess.Popen(["explorer.exe", win_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        elif sys_name == "Darwin":
            subprocess.Popen(["open", str(folder_path)])
        else:
            subprocess.Popen(["xdg-open", str(folder_path)])

    def open_output_folder(self):
        out_dir = Path(self.var_output.get().strip())
        out_dir.mkdir(parents=True, exist_ok=True)
        self.open_folder_in_os(out_dir)


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    app = DeepVirApp()
    app.mainloop()
