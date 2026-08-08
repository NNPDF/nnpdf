from collections import defaultdict
import dataclasses
import re

from .inspire import BibliographyEntry, fetch_bibtex

LATEX_ESCAPE = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}
_ESCAPE_RE = re.compile(
    r"\$[^$]*\$"  # Finds math mode substrings
    r"|\w+\\[a-zA-Z]+(?:\{[^}]*\})?|"  # Find latex_like strings
    + r"|".join(re.escape(c) for c in LATEX_ESCAPE)  # Find strings to escape
)
_EXTRACT_ID_RE = re.compile(r"/literature/(\d+)")

# Caption templates for grouped tables
_CAPTION_TEMPLATES = {
    "experiment": "from the {value} experiment",
    "process": "matching the {value} process",
    "nnpdf31_process": "in the {value} category (NNPDF 3.1)",
}


@dataclasses.dataclass
class LatexDatasetRow:
    """All information needed to write a LaTeX row."""

    dataset_name: str
    dataset_label: str
    observable_description: str
    experiment: str
    process: str
    nnpdf31_process: str
    bib_entry: BibliographyEntry

    @property
    def bib_key(self):
        return self.bib_entry.key

    @property
    def bibtex(self):
        return self.bib_entry.raw_bibtex

    def __post_init__(self):
        """Pass some strings thorugh escape-latex."""
        self.dataset_label = _escape_latex(self.dataset_label)
        self.observable_description = _escape_latex(self.observable_description)


def _escape_latex(text: str) -> str:
    """Escape LaTeX characters, unless they are completely within latex already.
    In that case they are matched (but not escaped).
    """  # https://xkcd.com/208/

    def skap(match):
        token = match.group()
        return LATEX_ESCAPE.get(token, token)

    return _ESCAPE_RE.sub(skap, text)


def build_latex_rows(datasets, group_by=None) -> dict:
    """Go over a list of datasets, use the inspire API to
    get their bibtex references and generate a list of LatexDatasetRow.

    If group_by is given, the output is grouped by the given key.
    Otherwise, the dictionary contain just one key: unsorted
    """
    rows = defaultdict(list)
    group_key = "unsorted"
    for dataset_name, dataset in datasets.items():
        inspire_url = dataset.inspire_url
        inspire_id = _EXTRACT_ID_RE.search(inspire_url).group(1)
        bib_entry = fetch_bibtex(inspire_id)
        row = LatexDatasetRow(
            dataset_name=dataset_name,
            dataset_label=dataset.plotting.dataset_label,
            observable_description=dataset.description,
            experiment=dataset.experiment,
            process=dataset.process,
            nnpdf31_process=dataset.nnpdf31_process,
            bib_entry=bib_entry,
        )
        if group_by is not None:
            group_key = getattr(dataset, group_by, "unsorted")
        rows[group_key].append(row)
    return rows


def generate_table(rows, group_by=None, group="unsorted"):
    """Given a list of rows, generate a LaTeX table.
    Generates a table with 3 columns:
        Dataset & Observable & Reference
    """
    if group_by is None:
        caption = "List of all datasets present in the runcard."
    else:
        caption = f"List of datasets matching {group_by}={group}"

    row_lines = []
    bibtex_list = []
    for row in rows:
        cite = f"\\cite{{{row.bib_key}}}"
        bibtex_list.append(row.bibtex)
        row_lines.append(" & ".join([row.dataset_label, row.observable_description, cite]))

    table_text = (
        "\n".join(
            [
                "\\begin{table}",
                "  \\centering",
                "  \\begin{tabular}{lll}",
                "    \\toprule",
                "    Dataset & Observable & Reference \\\\",
                "    \\midrule",
                *row_lines,
                "    \\bottomrule",
                "  \\end{tabular}",
                f"  \\caption{{{caption}}}",
                "  \\label{tab:datasets}",
                "\\end{table}",
            ]
        )
        + "\n"
    )
    return table_text, bibtex_list
