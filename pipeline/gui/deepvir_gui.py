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
import subprocess
import threading
import queue
import time
import webbrowser
import platform
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, font as tkfont

# Regular expression to match ANSI escape sequences (colors and cursor control)
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]|\r")

# =============================================================================
# CONSTANTS & GLOBAL PATH CONFIGURATION
# =============================================================================
APP_TITLE = "DeepVir - Viral Metagenomics & Dark Matter Discovery"
APP_VERSION = "v1.0 Desktop"

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
        "#                         version: 1.0  |  09.2026                                #",
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

        lbl_title = tk.Label(title_box, text="DeepVir", fg=COLOR_PRIMARY, bg=COLOR_CARD)
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

        lbl_dbs_header = tk.Label(card_dbs, text="External Databases (Optional)", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl_dbs_header, "head", "bold")
        lbl_dbs_header.pack(anchor="w", padx=14, pady=(12, 8))

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
        self.mod_reads_kraken = tk.BooleanVar(value=False)
        self.mod_reads_diamond = tk.BooleanVar(value=False)
        self.mod_assembly = tk.BooleanVar(value=True)
        self.mod_kraken2 = tk.BooleanVar(value=True)
        self.mod_diamond = tk.BooleanVar(value=True)
        self.mod_darkmatter = tk.BooleanVar(value=True)
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

        self.create_check_group(
            card_modules,
            "> Viral Discovery / Dark Matter Module",
            [
                ("RdRp catalytic motif detection via PalmAnnot (--darkmatter)", self.mod_darkmatter),
                ("Cross-validation against Reference Viral DB (--rvdb)", self.mod_rvdb),
            ],
        )

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
    # TAB 4: RESULTS & REPORTS BUILDER
    # =========================================================================
    def build_tab_results(self, parent):
        """
        Builds the contents of the Results & Reports tab (100% pure ASCII labels).
        """
        container = ttk.Frame(parent)
        container.pack(fill="both", expand=True, padx=6, pady=6)

        card = tk.Frame(container, bg=COLOR_CARD, highlightbackground=COLOR_CARD_BORDER, highlightthickness=1)
        card.pack(fill="both", expand=True, padx=12, pady=12)

        lbl = tk.Label(card, text="Results & Output Reports Viewer", fg=COLOR_ACCENT, bg=COLOR_CARD)
        self.register_scalable(lbl, "title", "bold")
        lbl.pack(anchor="w", padx=20, pady=(20, 6))

        lbl_sub = tk.Label(
            card,
            text="Open generated interactive reports directly in your default web browser:",
            fg=COLOR_MUTED,
            bg=COLOR_CARD,
        )
        self.register_scalable(lbl_sub, "body")
        lbl_sub.pack(anchor="w", padx=20, pady=(0, 20))

        # Action Buttons container
        btn_grid = tk.Frame(card, bg=COLOR_CARD)
        btn_grid.pack(fill="x", padx=20)

        # Open MultiQC report
        btn_mqc = ttk.Button(
            btn_grid,
            text="> Open MultiQC Quality Report (multiqc_report.html)",
            style="Secondary.TButton",
            command=lambda: self.open_output_file("multiqc_report.html"),
        )
        btn_mqc.pack(fill="x", pady=6)

        # Open Krona taxonomic wheels
        btn_krona = ttk.Button(
            btn_grid,
            text="> Open Krona Taxonomic Wheels (*krona*.html)",
            style="Secondary.TButton",
            command=lambda: self.open_output_pattern("*krona*.html"),
        )
        btn_krona.pack(fill="x", pady=6)

        # Open Dark Matter / Viral Diversity Report
        btn_dm = ttk.Button(
            btn_grid,
            text="> Open Viral Diversity / Dark Matter Report (*Report_Diversity*.html)",
            style="Secondary.TButton",
            command=lambda: self.open_output_pattern("*Report_Diversity*.html"),
        )
        btn_dm.pack(fill="x", pady=6)

        # Open Snakemake workflow summary report
        btn_smk = ttk.Button(
            btn_grid,
            text="> Open Snakemake Workflow Report (report.html)",
            style="Secondary.TButton",
            command=lambda: self.open_output_file("report.html"),
        )
        btn_smk.pack(fill="x", pady=6)

        # Open Centralized Execution Logs Folder (logs/)
        btn_logs = ttk.Button(
            btn_grid,
            text="> Open Execution Logs Folder (logs/)",
            style="Secondary.TButton",
            command=lambda: self.open_folder_in_os(LOGS_DIR),
        )
        btn_logs.pack(fill="x", pady=6)

        # Open Documentation & Licenses Folder (docs/)
        btn_docs = ttk.Button(
            btn_grid,
            text="> Open Documentation & Licenses Folder (docs/)",
            style="Secondary.TButton",
            command=lambda: self.open_folder_in_os(DOCS_DIR),
        )
        btn_docs.pack(fill="x", pady=6)

        # Open Output Folder in OS File Explorer
        btn_folder = ttk.Button(
            btn_grid,
            text="[OPEN] Output Folder in File Explorer (results/)",
            style="Primary.TButton",
            command=self.open_output_folder,
        )
        btn_folder.pack(fill="x", pady=16)

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
    # OUTPUT LAUNCHERS & SYSTEM INTEGRATION
    # =========================================================================
    def open_output_file(self, filename):
        out_dir = Path(self.var_output.get().strip())
        matches = list(out_dir.glob(f"**/{filename}"))
        if matches:
            target = matches[0]
            webbrowser.open(target.as_uri())
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
            webbrowser.open(matches[0].as_uri())
        else:
            messagebox.showwarning("Not Found", f"No files matching pattern '{pattern}' were found.")

    def open_folder_in_os(self, folder_path: Path):
        folder_path.mkdir(parents=True, exist_ok=True)
        sys_name = platform.system()
        if sys_name == "Windows":
            os.startfile(str(folder_path))
        elif sys_name == "Darwin":
            subprocess.Popen(["open", str(folder_path)])
        else:
            subprocess.Popen(["xdg-open", str(folder_path)])

    def open_output_folder(self):
        out_dir = Path(self.var_output.get().strip())
        if out_dir.exists():
            self.open_folder_in_os(out_dir)
        else:
            messagebox.showwarning("Warning", "The output directory does not exist yet.")


# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    app = DeepVirApp()
    app.mainloop()
