"""
blast_json_hsp_viewer_full.py

Full standalone script (Option A: stub functions).
Features:
 - Loads BLAST JSON (local file by default)
 - High-contrast PyQt5 GUI
 - Left: hits summary table
 - Right: HSP viewer with per-residue color highlighting (match=green, mismatch=red)
 - Slider to view windows of length 40 (adjustable)
 - Buttons: Isoform search, PTM search, Binding Info, General PTMs, Entire UniProt Page
   Each button opens a dynamic table window with results returned from test stub functions
 - All code heavily commented for clarity
"""

import os
import re
import json
from typing import Dict, Any, List, Optional

import requests  # used for URL-based JSON load option (safe even if unused)
import pandas as pd


from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTableWidget,
    QTableWidgetItem, QTextEdit, QSlider, QGroupBox, QScrollArea, QPushButton,
    QFileDialog, QMessageBox, QSizePolicy, QFrame, QTabWidget,
    QProgressDialog
)
from PyQt5.QtCore import Qt, QTimer

from pyteomics import mass

# -------------------------- Utils -----------------------------------
from .uniprot_querry import getFasta_from_accession, get_isoforms_list, get_isoforms_list, get_PTMs_from_accession_list, import_BindingInfo_uniProt, do_seq_Blast_get_json
from .general_utils import mw, enumerate_terminal_ptm_masses, fasta_to_sequence
# ---------------------------
# DEFAULT PATH (uploaded file)
# ---------------------------
#DEFAULT_JSON_PATH = "data/ncbiblast-R20251121-121303-0264-85327113-p2m.out.json"
DEFAULT_JSON_PATH = os.path.join(os.path.dirname(__file__), "data", "ncbiblast-R20251121-121303-0264-85327113-p2m.out.json")
# -----------------------------------------------------------------------------
# Helper: load JSON from local path or URL
# -----------------------------------------------------------------------------
def load_blast_json(path: Optional[str] = None, url: Optional[str] = None) -> Dict[str, Any]:
    """
    Load BLAST JSON either from a local file (path) or from a URL (url).
    - If both provided, URL takes precedence.
    - Returns parsed JSON dictionary.
    """
    if url:
        # Download JSON from URL
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        return r.json()

    if path:
        if not os.path.exists(path):
            raise FileNotFoundError(f"JSON file not found: {path}")
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    raise ValueError("Provide either path or url to load BLAST JSON.")

# -----------------------------------------------------------------------------
# Utility: replace spaces with '-' for alignment visibility
# -----------------------------------------------------------------------------
def replace_space_with_dash(s: Optional[str]) -> str:
    if s is None:
        return ""
    return s.replace(" ", "-")


def test_function_returning_json_url(user_input: str) -> str:
    """
    Dummy test function.
    You will replace this with your real BLAST request code.

    It receives the text from the top-panel text box
    and returns a JSON URL that the viewer will then load.
    """
    print("Test function called with:", user_input)

    json_url = do_seq_Blast_get_json(user_input)
    return json_url
    # simulate a returned JSON result URL
    #return("https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/ncbiblast-R20251121-121303-0264-85327113-p2m/json")

#
#
#

class BlastMainWindow(QWidget):
    """
    Wraps:
      - Top panel (JSON input + button)
      - BlastJsonHspViewer (bottom area)
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BLAST Viewer (with top panel)")
        self.resize(1500, 1000)

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ---------------- TOP PANEL ----------------
        top_row = QHBoxLayout()
        lbl = QLabel("Input:")
        lbl.setStyleSheet("color:white; font-weight:bold;")
        self.txt_input = QTextEdit()
        self.txt_input.setFixedHeight(40)
        self.txt_input.setPlaceholderText("Enter accession number or sequence here...")

        btn = QPushButton("Blast it wait a bit")
        btn.setStyleSheet("background:#444;color:white; font-weight:bold;")
        btn.clicked.connect(self._run_test_function)

        top_row.addWidget(lbl)
        top_row.addWidget(self.txt_input)
        top_row.addWidget(btn)

        layout.addLayout(top_row)

        # ---------------- HOLDER FOR EXISTING VIEWER ----------------
        self.viewer_container = QVBoxLayout()
        layout.addLayout(self.viewer_container)

        # Initially load default JSON
        self._load_new_viewer(DEFAULT_JSON_PATH, None)

    def _load_new_viewer(self, json_path=None, json_url=None):
        """
        Clears old viewer and loads a new one inside the bottom panel.
        """
        # remove old viewer (if exists)
        while self.viewer_container.count():
            item = self.viewer_container.takeAt(0)
            w = item.widget()
            if w:
                w.setParent(None)
                w.deleteLater()

        # load JSON
        blast_json = load_blast_json(path=json_path, url=json_url)

        # create viewer widget
        self.viewer = BlastJsonHspViewer(blast_json)
        self.viewer_container.addWidget(self.viewer)


    def _run_test_function(self):
        """
        Called when user presses the top-panel button.
        Shows a waiting dialog while test_function_returning_json_url()
        generates a JSON URL.
        """
        text = self.txt_input.toPlainText().strip()
        if not text:
            QMessageBox.warning(self, "Input needed", "Enter text before running test.")
            return

        # --------------------------
        # SHOW WAITING DIALOG
        # --------------------------
        wait = QProgressDialog("Running BLAST request...\nPlease wait.", None, 0, 0, self)
        wait.setWindowTitle("Please Wait")
        wait.setWindowModality(Qt.ApplicationModal)
        wait.setStyleSheet("background:#333;color:white;font-size:14px;")
        wait.setAutoClose(False)
        wait.setAutoReset(False)
        wait.show()

        # ---------------------------------------------
        # DELAYED EXECUTION (prevents GUI freezing)
        # ---------------------------------------------
        def run_and_update():
            try:
                # Call your test function (slow operation)
                json_url = test_function_returning_json_url(text)

                # Load new viewer
                self._load_new_viewer(json_url=json_url)

            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error running test function:\n{e}")

            # Close waiting dialog
            wait.close()

        # Process after 10 ms to allow dialog to show
        QTimer.singleShot(10, run_and_update)


# -----------------------------------------------------------------------------
# Re-usable dynamic table window (for button results)
# -----------------------------------------------------------------------------

class TableWindow(QWidget):
    """
    Dynamic table viewer.

    Supports both:
    - List of dictionaries  → column names from keys
    - List of lists         → optional headers; autogenerated if missing

    If list-of-lists is provided without headers, auto-create:
        ["Col1", "Col2", "Col3", ...]
    """

    def __init__(self, title: str, data: List, headers: Optional[List[str]] = None):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(800, 400)

        layout = QVBoxLayout()
        self.setLayout(layout)

        # Handle empty input
        if not data:
            layout.addWidget(QLabel("No data to display."))
            return

        # Determine row format
        first_row = data[0]
        is_dict = isinstance(first_row, dict)

        # -------------------------------
        # Case 1: LIST OF DICTIONARIES
        # -------------------------------
        if is_dict:
            columns = list(first_row.keys())
            rows = data

        # -------------------------------
        # Case 2: LIST OF LISTS
        # -------------------------------
        else:
            # Auto-create headers if none provided
            if headers is None:
                num_cols = len(first_row)
                headers = [f"Col{i+1}" for i in range(num_cols)]

            columns = headers
            rows = data

        # Build QTableWidget
        table = QTableWidget()
        table.setRowCount(len(rows))
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.setStyleSheet("background:#111;color:white;")

        # Fill data
        for r, row in enumerate(rows):
            if is_dict:
                # Dictionary-format rows
                for c, col in enumerate(columns):
                    value = row.get(col, "")
                    item = QTableWidgetItem(str(value))
                    item.setForeground(Qt.white)
                    table.setItem(r, c, item)

            else:
                # List-of-lists rows
                for c, value in enumerate(row):
                    item = QTableWidgetItem(str(value))
                    item.setForeground(Qt.white)
                    table.setItem(r, c, item)

        table.resizeColumnsToContents()
        layout.addWidget(table)

# -----------------------------------------------------------------------------
# USER CALLBACK FUNCTIONS
# -----------------------------------------------------------------------------
def user_query_info(qseq: str) -> str:
    """
    User-defined function that receives the full query sequence
    and returns a useful string to display beside the 'Query:' label.
    Replace this stub with your real function.
    """
    qseq_ = re.sub(r'[^A-Za-z]', '', qseq)

    #th_mass  = calculate_peptide_masses(qseq_)

    th_mass_mono = mw(qseq_, monoisotopic=True)
    th_mass_avg  = mw(qseq_, monoisotopic=False)
    #return f"Avg Mass={th_mass.get('average_mass')}  Monoisotopic Mass={th_mass.get('monoisotopic_mass')}"
    return f"Avg Mass={th_mass_avg}  Monoisotopic Mass={th_mass_mono}"

def user_subject_info(hseq: str) -> str:
    """
    User-defined function that receives the subject sequence
    and returns a useful string to display beside the 'Subject:' label.
    Replace this with your real function.
    """
    hseq_ = re.sub(r'[^A-Za-z]', '', hseq)    

    #th_mass  = calculate_peptide_masses(hseq_)
    #return f"Avg Mass={th_mass.get('average_mass')}  Monoisotopic Mass={th_mass.get('monoisotopic_mass')}"

    th_mass_mono = mw(hseq_, monoisotopic=True)
    th_mass_avg  = mw(hseq_, monoisotopic=False)
    return f"Avg Mass={th_mass_avg}  Monoisotopic Mass={th_mass_mono}"


# -----------------------------------------------------------------------------
# HSP viewer widget (MODIFIED)
# -----------------------------------------------------------------------------
class HspViewerWidget(QWidget):
    """
    Widget that displays HSP metadata and a sliding, colorized view (window)
    of hsp_qseq, hsp_mseq, hsp_hseq. Matches=green, mismatches=red.
    Modified to show user-function output next to Query/Subject labels.
    """

    def __init__(self, hsp: Dict[str, Any], window_size: int = 40, parent=None):
        super().__init__(parent)
        self.hsp = hsp
        self.window_size = window_size

        # prepare sequences
        self.qseq = replace_space_with_dash(hsp.get("hsp_qseq", "") or "")
        self.mseq = replace_space_with_dash(hsp.get("hsp_mseq", "") or "")
        self.hseq = replace_space_with_dash(hsp.get("hsp_hseq", "") or "")

        # pad to same length
        max_len = max(len(self.qseq), len(self.mseq), len(self.hseq))
        self.qseq = self.qseq.ljust(max_len, "-")
        self.mseq = self.mseq.ljust(max_len, "-")
        self.hseq = self.hseq.ljust(max_len, "-")
        self.seq_length = max_len

        # compute user-supplied info strings
        self.query_info_text = user_query_info(self.qseq)
        self.subject_info_text = user_subject_info(self.hseq)

        self._build_ui()
        self._update_display(0)

    def _build_ui(self):
        v = QVBoxLayout()
        v.setSpacing(6)

        # metadata row
        meta_row = QHBoxLayout()
        for key, label in [
            ("hsp_num", "HSP #"), ("hsp_bit_score", "Bit score"),
            ("hsp_score", "Score"), ("hsp_expect", "E-value"),
            ("hsp_identity", "Identity"), ("hsp_align_len", "Align len"),
            ("hsp_gaps", "Gaps")
        ]:
            val = self.hsp.get(key, "")
            lbl = QLabel(f"<b>{label}:</b> {val}")
            lbl.setStyleSheet("color:white;")
            meta_row.addWidget(lbl)
        meta_row.addStretch()
        v.addLayout(meta_row)

        # text style
        text_style = (
            "font-family:monospace; font-size:11pt; "
            "background:#111; color:white;"
        )

        # ----- Query label + user info -----
        qrow = QHBoxLayout()
        qlabel = QLabel(f"<span style='color:white'><b>Query:</b></span>")
        qinfo = QLabel(f"<span style='color:lightblue'>{self.query_info_text}</span>")
        qrow.addWidget(qlabel)
        qrow.addWidget(qinfo)
        qrow.addStretch()
        v.addLayout(qrow)

        self.qtext = QTextEdit()
        self.qtext.setReadOnly(True)
        self.qtext.setLineWrapMode(QTextEdit.NoWrap)
        self.qtext.setFixedHeight(34)
        self.qtext.setStyleSheet(text_style)
        v.addWidget(self.qtext)

        # ----- Match -----
        v.addWidget(QLabel("<span style='color:white'><b>Match:</b></span>"))
        self.mtext = QTextEdit()
        self.mtext.setReadOnly(True)
        self.mtext.setLineWrapMode(QTextEdit.NoWrap)
        self.mtext.setFixedHeight(34)
        self.mtext.setStyleSheet(text_style)
        v.addWidget(self.mtext)

        # ----- Subject label + user info -----
        srow = QHBoxLayout()
        slabel = QLabel("<span style='color:white'><b>Subject:</b></span>")
        sinfo = QLabel(f"<span style='color:lightblue'>{self.subject_info_text}</span>")
        srow.addWidget(slabel)
        srow.addWidget(sinfo)
        srow.addStretch()
        v.addLayout(srow)

        self.htext = QTextEdit()
        self.htext.setReadOnly(True)
        self.htext.setLineWrapMode(QTextEdit.NoWrap)
        self.htext.setFixedHeight(34)
        self.htext.setStyleSheet(text_style)
        v.addWidget(self.htext)

        # Slider + position
        slider_row = QHBoxLayout()
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(max(0, self.seq_length - self.window_size))
        self.slider.valueChanged.connect(self._update_display)
        slider_row.addWidget(self.slider)

        self.pos_label = QLabel("")
        self.pos_label.setStyleSheet("color:white;")
        slider_row.addWidget(self.pos_label)
        v.addLayout(slider_row)

        self.setLayout(v)
        self.setStyleSheet("background:#222;")

    def _colorize_html(self, qwin, mwin, hwin):
        q_html, m_html, h_html = [], [], []
        for qc, mc, hc in zip(qwin, mwin, hwin):
            color = "lightgreen" if qc == hc else "#ff5555"
            q_html.append(f"<span style='color:{color}'>{qc}</span>")
            m_html.append(f"<span style='color:{color}'>{mc}</span>")
            h_html.append(f"<span style='color:{color}'>{hc}</span>")
        return "".join(q_html), "".join(m_html), "".join(h_html)

    def _update_display(self, pos_value):
        pos = max(0, min(int(pos_value), max(0, self.seq_length - self.window_size)))
        end = pos + self.window_size

        q_html, m_html, h_html = self._colorize_html(
            self.qseq[pos:end], self.mseq[pos:end], self.hseq[pos:end]
        )

        self.qtext.setHtml(q_html)
        self.mtext.setHtml(m_html)
        self.htext.setHtml(h_html)
        self.pos_label.setText(
            f"<span style='color:white'>{pos+1}-{end} of {self.seq_length}</span>"
        )

# -----------------------------------------------------------------------------
# Test stub functions for button actions (Option A)
# Replace these with real API calls in future if desired
# -----------------------------------------------------------------------------
def test_isoform_search(hit_accession: str) -> List[Dict[str, Any]]:
    """Return demo isoform rows for a given accession."""
    return get_isoforms_list(hit_accession)
    #return [
    #    {"Isoform": f"{hit_accession}-1", "Length": 410, "Evidence": "Experimental"},
    #    {"Isoform": f"{hit_accession}-2", "Length": 387, "Evidence": "Predicted"},
    #]

def test_ptm_search(hit_accession: str) -> List[Dict[str, Any]]:
    """Return demo PTM rows."""
    return get_PTMs_from_accession_list(hit_accession)
    #return [
    #    {"Position": 53, "Type": "Phosphorylation", "Source": "UniProt"},
    #    {"Position": 129, "Type": "Acetylation", "Source": "UniProt"},
    #]

def test_binding_info(hit_accession: str) -> List[Dict[str, Any]]:
    """Return demo binding site rows."""
    return import_BindingInfo_uniProt(hit_accession)
    #return [
    #    {"Ligand": "ATP", "Position": 201, "Evidence": "Experimental"},
    #    {"Ligand": "Zn2+", "Position": 349, "Evidence": "Modeled"},
    #]

def test_general_info(hit_accession: str) -> List[Dict[str, Any]]:
    """Return demo general info rows."""
    fasta_seq = getFasta_from_accession(hit_accession)
    seq       = fasta_to_sequence(fasta_seq)
    return(enumerate_terminal_ptm_masses(seq))
#    return [
#        {"Field": "Protein Name", "Value": "Example Kinase"},
#        {"Field": "Organism", "Value": "Homo sapiens"},
#        {"Field": "Mass (Da)", "Value": "53,450"},
#    ]

def test_entire_uniprot_page(hit_accession: str) -> List[Dict[str, Any]]:
    """Return demo extracted UniProt page sections."""
    return [
        {"Section": "Function", "Text": "Catalyzes transfer of phosphate groups."},
        {"Section": "Pathway", "Text": "Glycolysis"},
        {"Section": "Subcellular location", "Text": "Cytoplasm"},
    ]


# -------------------------------------------------------------------------
# UniProt JSON fetch helper
# -------------------------------------------------------------------------
def fetch_uniprot_json(accession: str) -> dict:
    """
    Query the UniProt REST API to fetch a full JSON record for an accession.

    Parameters
    ----------
    accession : str
        UniProt accession (e.g. 'P31749').

    Returns
    -------
    dict
        Parsed UniProt JSON.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}?format=json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


# -------------------------------------------------------------------------
# Full UniProt browser window (tabbed)
# -------------------------------------------------------------------------
class UniProtBrowserWindow(QWidget):
    """
    A richer UniProt "browser" for a single accession.

    Tabs:
      - Overview
      - Sequence
      - Features
      - PTMs
      - GO / Ontology
      - Cross-refs
    """

    def __init__(self, accession: str, parent=None):
        super().__init__(parent)
        self.accession = accession
        self.setWindowTitle(f"UniProt Browser: {accession}")
        self.resize(1000, 750)

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(8, 8, 8, 8)

        # Top title label
        self.title_label = QLabel(f"<h2 style='color:#00e0ff'>UniProt: {accession}</h2>")
        main_layout.addWidget(self.title_label)

        # Tab widget
        self.tabs = QTabWidget()
        self.tabs.setStyleSheet("QTabWidget::pane { border: 1px solid #444; } "
                                "QTabBar::tab { background:#222; color:white; padding:4px 10px;} "
                                "QTabBar::tab:selected { background:#444; }")
        main_layout.addWidget(self.tabs)

        # Fetch data and populate tabs
        try:
            self.data = fetch_uniprot_json(accession)
        except Exception as e:
            err = QTextEdit()
            err.setReadOnly(True)
            err.setStyleSheet("background:#111;color:#ff5555;")
            err.setText(f"Error fetching UniProt entry {accession}:\n{e}")
            self.tabs.addTab(err, "Error")
            return

        # Create all tabs
        self._build_overview_tab()
        self._build_sequence_tab()
        self._build_features_tab()
        self._build_ptm_tab()
        self._build_go_tab()
        self._build_crossref_tab()

        self.setStyleSheet("background:#111;")

    # ------------------------- helpers -------------------------
    def _safe(self, d, path, default=None):
        cur = d
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                return default
            cur = cur[k]
        return cur

    # ------------------------- tab builders -------------------------
    def _build_overview_tab(self):
        """
        Overview: protein name, organism, gene names, length, function, etc.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        text = QTextEdit()
        text.setReadOnly(True)
        text.setLineWrapMode(QTextEdit.WidgetWidth)
        text.setStyleSheet("background:#111;color:white;font-family:monospace;font-size:11pt;")

        # Protein name
        pname = self._safe(self.data, ["proteinDescription", "recommendedName", "fullName", "value"], "")
        # Organism
        org = self._safe(self.data, ["organism", "scientificName"], "")
        # Length
        length = self._safe(self.data, ["sequence", "length"], "")
        # Compute mass if possible
        seq = self._safe(self.data, ["sequence", "value"], "")
        mass_val = None
        if seq:
            try:
                mass_val = mass.calculate_mass(seq)
            except Exception:
                mass_val = None

        # Gene names
        gene_names = []
        for g in self.data.get("genes", []):
            if "geneName" in g:
                gene_names.append(g["geneName"]["value"])
            for syn in g.get("synonyms", []):
                gene_names.append(syn["value"])
        gene_names_str = ", ".join(sorted(set(gene_names))) if gene_names else "N/A"

        # Function comment(s)
        func_texts = []
        for c in self.data.get("comments", []):
            if c.get("type") == "FUNCTION":
                for t in c.get("texts", []):
                    func_texts.append(t.get("value", ""))
        func_html = "<br>".join(func_texts) if func_texts else "N/A"

        overview_html = []
        overview_html.append("<div style='color:white'>")
        overview_html.append(f"<p><b style='color:#ffd700'>Accession:</b> {self.accession}</p>")
        overview_html.append(f"<p><b style='color:#ffd700'>Protein:</b> {pname}</p>")
        overview_html.append(f"<p><b style='color:#ffd700'>Organism:</b> {org}</p>")
        overview_html.append(f"<p><b style='color:#ffd700'>Gene names:</b> {gene_names_str}</p>")
        overview_html.append(f"<p><b style='color:#ffd700'>Length:</b> {length}</p>")
        if mass_val is not None:
            overview_html.append(f"<p><b style='color:#ffd700'>Theoretical mass:</b> {mass_val:.2f} Da</p>")

        overview_html.append("<hr style='border:1px solid #444;'>")
        overview_html.append("<p><b style='color:#00e0ff'>Function:</b><br>")
        overview_html.append(f"<span style='color:#e0e0e0'>{func_html}</span></p>")
        overview_html.append("</div>")

        text.setHtml("\n".join(overview_html))
        v.addWidget(text)
        self.tabs.addTab(w, "Overview")

    def _build_sequence_tab(self):
        """
        Sequence tab: raw sequence with numbering, monospaced, scrollable in X/Y.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        seq = self._safe(self.data, ["sequence", "value"], "")
        if not seq:
            txt = QTextEdit("No sequence available.")
            txt.setReadOnly(True)
            txt.setStyleSheet("background:#111;color:white;")
            v.addWidget(txt)
            self.tabs.addTab(w, "Sequence")
            return

        # Format sequence with numbering every 10 residues
        lines = []
        block = 60
        for i in range(0, len(seq), block):
            chunk = seq[i:i+block]
            lines.append(f"{i+1:6d}  {chunk}")

        txt = QTextEdit()
        txt.setReadOnly(True)
        txt.setLineWrapMode(QTextEdit.NoWrap)  # horizontal scroll
        txt.setStyleSheet("background:#111;color:#e0e0e0;font-family:monospace;font-size:11pt;")
        txt.setText("\n".join(lines))

        v.addWidget(txt)
        self.tabs.addTab(w, "Sequence")

    def _build_features_tab(self):
        """
        Features tab: UniProt 'features' array: regions, domains, sites, etc.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        feats = self.data.get("features", []) or []
        table = QTableWidget()
        table.setStyleSheet("background:#111;color:white;gridline-color:#444;")
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Type", "Description", "Begin", "End", "Ligand"])

        table.setRowCount(len(feats))
        for r, ft in enumerate(feats):
            ftype = ft.get("type", "")
            desc = ft.get("description", "")
            loc = ft.get("location", {})
            beg = loc.get("start", {}).get("value", "")
            end = loc.get("end", {}).get("value", beg)
            ligand = ft.get("ligand", {}).get("name", "")

            for c, val in enumerate([ftype, desc, str(beg), str(end), ligand]):
                item = QTableWidgetItem(str(val))
                item.setForeground(Qt.white)
                table.setItem(r, c, item)

        table.resizeColumnsToContents()
        table.setHorizontalScrollMode(QTableWidget.ScrollPerPixel)
        table.setVerticalScrollMode(QTableWidget.ScrollPerPixel)
        v.addWidget(table)
        self.tabs.addTab(w, "Features")

    def _build_ptm_tab(self):
        """
        PTMs tab: combine PTM-related comments and modified residue features.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        rows = []

        # PTM comments
        for c in self.data.get("comments", []):
            if c.get("type") in ("PTM", "MODIFIED RESIDUE"):
                for t in c.get("texts", []):
                    rows.append({
                        "Source": "Comment",
                        "Position": "",
                        "Type": c.get("type", ""),
                        "Detail": t.get("value", "")
                    })

        # Modified residue features
        for ft in self.data.get("features", []):
            if ft.get("type") in ("MODIFIED RESIDUE", "LIPIDATION", "GLYCOSYLATION"):
                loc = ft.get("location", {})
                pos = loc.get("start", {}).get("value", "")
                desc = ft.get("description", "")
                rows.append({
                    "Source": "Feature",
                    "Position": pos,
                    "Type": ft.get("type", ""),
                    "Detail": desc
                })

        if not rows:
            txt = QTextEdit("No PTM information available.")
            txt.setReadOnly(True)
            txt.setStyleSheet("background:#111;color:white;")
            v.addWidget(txt)
            self.tabs.addTab(w, "PTMs")
            return

        table = QTableWidget()
        table.setStyleSheet("background:#111;color:white;gridline-color:#444;")
        cols = ["Source", "Position", "Type", "Detail"]
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.setRowCount(len(rows))

        for r, row in enumerate(rows):
            for c, col in enumerate(cols):
                item = QTableWidgetItem(str(row.get(col, "")))
                item.setForeground(Qt.white)
                table.setItem(r, c, item)

        table.resizeColumnsToContents()
        v.addWidget(table)
        self.tabs.addTab(w, "PTMs")

    def _build_go_tab(self):
        """
        GO / Ontology tab: GO ID, term, aspect.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        rows = []
        for xref in self.data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "GO":
                go_id = xref.get("id", "")
                term = ""
                aspect = ""
                for p in xref.get("properties", []):
                    if p["key"] == "term":
                        term = p["value"]
                    if p["key"] == "aspect":
                        aspect = p["value"]
                rows.append({"GO ID": go_id, "Aspect": aspect, "Term": term})

        if not rows:
            txt = QTextEdit("No GO / ontology information available.")
            txt.setReadOnly(True)
            txt.setStyleSheet("background:#111;color:white;")
            v.addWidget(txt)
            self.tabs.addTab(w, "GO / Ontology")
            return

        table = QTableWidget()
        table.setStyleSheet("background:#111;color:white;gridline-color:#444;")
        cols = ["GO ID", "Aspect", "Term"]
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.setRowCount(len(rows))

        for r, row in enumerate(rows):
            for c, col in enumerate(cols):
                item = QTableWidgetItem(str(row.get(col, "")))
                item.setForeground(Qt.white)
                table.setItem(r, c, item)

        table.resizeColumnsToContents()
        v.addWidget(table)
        self.tabs.addTab(w, "GO / Ontology")

    def _build_crossref_tab(self):
        """
        Cross-refs tab: all UniProtKBCrossReferences - DB, ID, extra properties.
        """
        w = QWidget()
        v = QVBoxLayout(w)

        xrefs = self.data.get("uniProtKBCrossReferences", []) or []
        if not xrefs:
            txt = QTextEdit("No cross-reference information available.")
            txt.setReadOnly(True)
            txt.setStyleSheet("background:#111;color:white;")
            v.addWidget(txt)
            self.tabs.addTab(w, "Cross-refs")
            return

        table = QTableWidget()
        table.setStyleSheet("background:#111;color:white;gridline-color:#444;")
        cols = ["Database", "Accession", "Properties"]
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.setRowCount(len(xrefs))

        for r, xr in enumerate(xrefs):
            db = xr.get("database", "")
            acc = xr.get("id", "")
            props = "; ".join(
                f"{p.get('key')}={p.get('value')}"
                for p in xr.get("properties", [])
            )
            for c, val in enumerate([db, acc, props]):
                item = QTableWidgetItem(str(val))
                item.setForeground(Qt.white)
                table.setItem(r, c, item)

        table.resizeColumnsToContents()
        v.addWidget(table)
        self.tabs.addTab(w, "Cross-refs")


# -----------------------------------------------------------------------------
# Main GUI: hits table on left, HSP viewer + buttons on right
# -----------------------------------------------------------------------------
class BlastJsonHspViewer(QWidget):
    """
    Main application window: left = hits summary table, right = HSP viewer area.
    Buttons above HSPs open dynamic table windows with results from stub functions.
    """
    def __init__(self, blast_json: Dict[str, Any], window_size: int = 40):
        super().__init__()
        self.setWindowTitle("BLAST JSON HSP Viewer (Stub functions)")
        self.resize(1400, 900)
        self.window_size = window_size

        # parse hits from JSON; JSON expected to have top-level "hits" list
        # structure: {"hits": [ {hit fields... , "hit_hsps": [hsp dicts...] }, ... ]}
        self.data = blast_json
        self.hits = blast_json.get("hits", []) or []

        self._build_ui()

    def _build_ui(self):
        # outer layout
        root_layout = QHBoxLayout()
        root_layout.setContentsMargins(8, 8, 8, 8)
        root_layout.setSpacing(12)

        # ---------- Left panel: hits summary table ----------
        left_layout = QVBoxLayout()
        left_label = QLabel("<b style='color:white'>Hits (select row)</b>")
        left_label.setStyleSheet("color:white;")
        left_layout.addWidget(left_label)

        # choose columns to show in summary table (best-effort)
        cols = ["hit_num", "hit_acc", "hit_id", "hit_def", "hit_os", "hit_len"]
        table = QTableWidget()
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.setRowCount(len(self.hits))
        table.setStyleSheet("background:#111;color:white;gridline-color:#444;")
        table.verticalHeader().setVisible(False)

        for r, hit in enumerate(self.hits):
            for c, col in enumerate(cols):
                val = hit.get(col, "")
                item = QTableWidgetItem(str(val))
                item.setForeground(Qt.white)
                item.setBackground(Qt.black)
                table.setItem(r, c, item)
        table.resizeColumnsToContents()
        table.setSelectionBehavior(QTableWidget.SelectRows)
        table.cellClicked.connect(self._on_hit_selected)
        left_layout.addWidget(table)

        # export button
        btn_export = QPushButton("Export summary CSV")
        btn_export.setStyleSheet("background:#2b7;color:black;")
        btn_export.clicked.connect(self._export_summary_csv)
        left_layout.addWidget(btn_export)

        left_widget = QWidget()
        left_widget.setLayout(left_layout)
        left_widget.setFixedWidth(520)

        # ---------- Right panel: buttons + HSP scroll area ----------
        right_layout = QVBoxLayout()
        right_header = QLabel("<b style='color:white'>HSPs & Tools</b>")
        right_header.setStyleSheet("color:white;")
        right_layout.addWidget(right_header)

        # button row
        btn_row = QHBoxLayout()
        self.btn_isoform = QPushButton("Isoform search")
        self.btn_ptm = QPushButton("PTM search")
        self.btn_bind = QPushButton("Binding Info")
        self.btn_general = QPushButton("General PTMs")
        self.btn_uniprot = QPushButton("Entire uniprot page")
        for b in (self.btn_isoform, self.btn_ptm, self.btn_bind, self.btn_general, self.btn_uniprot):
            b.setStyleSheet("background:#333;color:white;font-weight:bold;")
            btn_row.addWidget(b)
        right_layout.addLayout(btn_row)

        # connect buttons to handlers
        self.btn_isoform.clicked.connect(self._do_isoform)
        self.btn_ptm.clicked.connect(self._do_ptm)
        self.btn_bind.clicked.connect(self._do_binding)
        self.btn_general.clicked.connect(self._do_general)
        self.btn_uniprot.clicked.connect(self._do_uniprot)

        # scroll area for HSP widgets
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.hsp_list_container = QWidget()
        self.hsp_list_layout = QVBoxLayout()
        self.hsp_list_layout.setSpacing(10)
        self.hsp_list_layout.addStretch()
        self.hsp_list_container.setLayout(self.hsp_list_layout)
        self.hsp_list_container.setStyleSheet("background:#222;")
        self.scroll_area.setWidget(self.hsp_list_container)

        right_layout.addWidget(self.scroll_area)
        right_widget = QWidget()
        right_widget.setLayout(right_layout)

        # assemble
        root_layout.addWidget(left_widget)
        root_layout.addWidget(right_widget, 1)
        self.setLayout(root_layout)
        self.table = table

        self.setStyleSheet("background:#111;")

        # keep references to open TableWindows (prevent GC)
        self._open_windows = []

    # -------------------------
    # Handlers for selecting a hit
    # -------------------------
    def _on_hit_selected(self, row: int, col: int):
        """
        When user clicks a row in hits table, populate the right scroll area
        with HSP widgets for that hit.
        """
        # clear existing HSP widgets
        while self.hsp_list_layout.count():
            item = self.hsp_list_layout.takeAt(0)
            if item is None:
                break
            w = item.widget()
            if w is not None:
                w.deleteLater()

        if row < 0 or row >= len(self.hits):
            return

        hit = self.hits[row]
        # header for the hit
        hdr = QLabel(f"<span style='color:#ffd700'><b>Hit {hit.get('hit_num', '')}:</b></span> "
                     f"<span style='color:white'>{hit.get('hit_acc', '')} - {hit.get('hit_def','')}</span>")
        hdr.setStyleSheet("background:#333;padding:6px;border-radius:4px;")
        self.hsp_list_layout.addWidget(hdr)

        hsps = hit.get("hit_hsps", []) or []
        if not hsps:
            no_lbl = QLabel("<span style='color:white'>No HSPs found for this hit.</span>")
            self.hsp_list_layout.addWidget(no_lbl)
        else:
            for hsp in hsps:
                box = QGroupBox(f"HSP #{hsp.get('hsp_num','')}")
                box.setStyleSheet("QGroupBox { color: white; border: 1px solid #555; margin-top: 6px; padding: 6px; }")
                inner = QVBoxLayout()
                hv = HspViewerWidget(hsp, window_size=self.window_size)
                inner.addWidget(hv)
                box.setLayout(inner)
                self.hsp_list_layout.addWidget(box)

        self.hsp_list_layout.addStretch()

    # -------------------------
    # Button actions: use test stubs and open TableWindow
    # -------------------------
    def _get_selected_hit_accession(self) -> Optional[str]:
        selected_items = self.table.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No selection", "Please select a hit row first.")
            return None
        row = selected_items[0].row()
        hit = self.hits[row]
        return hit.get("hit_acc")

    def _open_table_window(self, title: str, rows: List[Dict[str, Any]]):
        win = TableWindow(title, rows)
        win.show()
        # keep a reference
        self._open_windows.append(win)

    def _safe_call(self, title: str, func, accession: str):
        """
        Helper wrapper so that API exceptions never kill the GUI.
        Shows a QMessageBox with the error.
        """
        try:
           rows = func(accession)
           if not isinstance(rows, list):
                raise ValueError(f"{title} returned a non-list object.")
           self._open_table_window(title, rows)

        except Exception as e:
           QMessageBox.critical(
                self,
                f"{title} Error",
                f"An error occurred while performing '{title}':\n\n{str(e)}"
            )

    def _do_isoform(self):
        acc = self._get_selected_hit_accession()
        if not acc:
            return
        rows = test_isoform_search(acc)
        #self._open_table_window("Isoform Search", rows)
        self._safe_call("Isoform Search", test_isoform_search, acc)

    def _do_ptm(self):
        acc = self._get_selected_hit_accession()
        if not acc:
            return
        rows = test_ptm_search(acc)
        #self._open_table_window("PTM Search", rows)
        self._safe_call("PTM Search", test_ptm_search, acc)

    def _do_binding(self):
        acc = self._get_selected_hit_accession()
        if not acc:
            return
        rows = test_binding_info(acc)
        #self._open_table_window("Binding Info", rows)
        self._safe_call("Binding Info", test_binding_info, acc)

    def _do_general(self):
        acc = self._get_selected_hit_accession()
        if not acc:
            return
        rows = test_general_info(acc)
        #self._open_table_window("General PTMs", rows)
        self._safe_call("General PTMs", test_general_info, acc)

    #def _do_uniprot(self):
    #    acc = self._get_selected_hit_accession()
    #    if not acc:
    #        return
    #    rows = test_entire_uniprot_page(acc)
    #    self._open_table_window("Entire UniProt Page", rows)
    
    def _do_uniprot(self):
         """
         Open the full UniProt browser window for the selected accession.
         Shows a small wait dialog while fetching JSON.
         """
         acc = self._get_selected_hit_accession()
         if not acc:
            return

         wait = QProgressDialog(f"Fetching UniProt entry {acc}...", None, 0, 0, self)
         wait.setWindowTitle("Please wait")
         wait.setWindowModality(Qt.ApplicationModal)
         wait.setStyleSheet("background:#333;color:white;")
         wait.setAutoClose(False)
         wait.setAutoReset(False)
         wait.show()

         def open_browser():
            try:
               win = UniProtBrowserWindow(acc)
               win.show()
               # keep a reference alive
               self._open_windows.append(win)
            except Exception as e:
               QMessageBox.critical(self, "Error", f"Failed to open UniProt browser:\n{e}")
            wait.close()

         QTimer.singleShot(20, open_browser)

    
    # -------------------------
    # Export summary to CSV
    # -------------------------
    def _export_summary_csv(self):
        rows = []
        for hit in self.hits:
            rows.append({
                "hit_num": hit.get("hit_num"),
                "hit_acc": hit.get("hit_acc"),
                "hit_id": hit.get("hit_id"),
                "hit_def": hit.get("hit_def"),
                "hit_os": hit.get("hit_os"),
                "hit_len": hit.get("hit_len"),
                "top_hsp_bit": (hit.get("hit_hsps") or [{}])[0].get("hsp_bit_score") if hit.get("hit_hsps") else None,
                "top_hsp_e": (hit.get("hit_hsps") or [{}])[0].get("hsp_expect") if hit.get("hit_hsps") else None,
            })
        df = pd.DataFrame(rows)

        fname, _ = QFileDialog.getSaveFileName(self, "Save summary CSV", "blast_summary.csv", "CSV files (*.csv)")
        if not fname:
            return
        try:
            df.to_csv(fname, index=False)
            QMessageBox.information(self, "Saved", f"Summary saved to {fname}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save CSV: {e}")

# -----------------------------------------------------------------------------
# Top-level run function (public API)
# -----------------------------------------------------------------------------
def run_blast_hsp_viewer(json_path: Optional[str] = DEFAULT_JSON_PATH, json_url: Optional[str] = None, window_size: int = 40):
    """
    Load JSON (path or url) and run the GUI.
    Default uses the uploaded file path.
    """
    blast_data = load_blast_json(path=json_path, url=json_url)

    #app = QApplication([])
    # ❗ IMPORTANT FIX:
    app = QApplication.instance()
    created_app = False
    
    if app is None:
        app = QApplication([])
        created_app = True

    win = BlastMainWindow()
    win.show()

    #gui = BlastJsonHspViewer(blast_data, window_size=window_size)
    #gui.show()
    # Only start the loop if we created this QApplication
    if created_app:
        app.exec_()
    #app.exec_()

# -----------------------------------------------------------------------------
# Simple tests for stub functions
# -----------------------------------------------------------------------------
def test_stubs():
    print("Testing stub functions with accession 'P12345':")
    print("Isoform:", test_isoform_search("P12345"))
    print("PTMs:", test_ptm_search("P12345"))
    print("Binding:", test_binding_info("P12345"))
    print("General:", test_general_info("P12345"))
    print("UniProt page:", test_entire_uniprot_page("P12345"))

# -----------------------------------------------------------------------------
# If run as a script, launch GUI with default path; also run stub tests in console.
# -----------------------------------------------------------------------------
#if __name__ == "__main__":
def main():
    # Option A: local bundled JSON
    json_path = os.path.join(os.path.dirname(__file__), "data", "ncbiblast-R20251121-121303-0264-85327113-p2m.out.json")
    run_blast_hsp_viewer(json_path=json_path)


#    # Run stub tests in the console first
#    #test_stubs()
#
#    # Then open the GUI using the default uploaded JSON
#    try:
#        run_blast_hsp_viewer(json_url="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/ncbiblast-R20251121-121303-0264-85327113-p2m/json")
#    except FileNotFoundError as e:
#        print("ERROR:", e)
#        print("If the default JSON path is not present, supply a local path or URL to run_blast_hsp_viewer().")


if __name__ == "__main__":
    main()
