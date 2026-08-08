"""
Entry point for the nnpdf-data CLI with the following options:

nnpdf-data latex <runcard.yml>
    Takes a NNPDF runcard and generate all biblatex entry in its ``dataset_inputs`` field
    toghether with latex tables.

nnpdf-data list <filter> [--runcard]
    List of all datasets matching the glob fliter <filter>
"""

import argparse
import logging
from pathlib import Path

from . import index, latex
from ..utils import yaml_safe


def _setup_logging(verbose=False):
    """Configure logging based on verbosity."""
    level = logging.DEBUG if verbose else logging.WARNING
    fmt = "%(levelname)s: %(name)s: %(message)s" if verbose else "%(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=fmt)


def _print_table(entries) -> None:
    """Print list entries with aligned columns."""
    exp_w = max((len(e.experiment) for e in entries.values()), default=len("EXPERIMENT"))
    proc_w = max((len(e.process) for e in entries.values()), default=len("PROCESS"))

    print(f"{'EXPERIMENT':<{exp_w}}  {'PROCESS':<{proc_w}}  DATASET")
    for name, e in entries.items():
        print(f"{e.experiment:<{exp_w}}  {e.process:<{proc_w}}  {name}")


def _cmd_list(entries, filter_pattern=None, sort_mode=None, runcard=False):
    """Run the list command."""
    if filter_pattern is not None:
        entries = index.filter_index(entries, filter_pattern)
        if not entries:
            print(f"No entries found with filter: {filter_pattern}")
            return -1

    if sort_mode is not None:
        entries = index.sort_index(entries, sort_mode)

    if runcard:
        for entry in entries.keys():
            print(f"- {{ dataset: {entry} }}")
    else:
        _print_table(entries)


def _cmd_latex(entries, runcard_path, sort_mode=None, group_by=None, output_bib=None):
    """Run the latex command."""
    if not runcard_path.is_file():
        raise FileNotFoundError(f"Runcard not found at {runcard_path}")

    runcard_info = yaml_safe.load(runcard_path.read_text(encoding="utf-8"))
    dataset_inputs = runcard_info["dataset_inputs"]

    datasets = {}
    for dataset_input in dataset_inputs:
        dataset_name = dataset_input["dataset"]
        try:
            datasets[dataset_name] = entries[dataset_name]
        except KeyError as e:
            raise KeyError(
                f"{dataset_name} from the runcard is not found among the available dataset"
            ) from e
    # At this point datasets is ordered as they are in the runcard

    if sort_mode is not None:
        datasets = index.sort_index(datasets, sort_mode)

    group_by_key = None
    if group_by:
        group_by_key = sort_mode

    # Latex rows is a dictionary of {group: [rows]}
    latex_rows = latex.build_latex_rows(datasets, group_by=group_by)
    latex_printout = []
    bibtex_entries = []
    for group, row_list in latex_rows.items():
        table, bibtex = latex.generate_table(row_list, group_by=group_by, group=group)
        latex_printout.append(table)
        bibtex_entries += bibtex

    print("\n\n".join(latex_printout))

    bibtex_text = "\n".join(bibtex_entries)

    if output_bib is None:
        print(bibtex_text)
    else:
        output_bib.write_text("\n".join(bibtex_entries))


def main():
    parser = argparse.ArgumentParser(prog="nnpdf-data", description="NNPDF data utilities.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    list_parser = subparsers.add_parser("list", help="List all available datasets.")
    list_parser.add_argument(
        "filter",
        nargs="?",
        default="*",
        help="Glob pattern to filter dataset. Supports brace expansion (e.g., '{ATLAS,CMS}*')p",
    )
    list_parser.add_argument(
        "--yaml",
        action="store_true",
        help="Output in the NNPDF yaml format ready for ``dataset_inputs:``",
    )

    latex_parser = subparsers.add_parser(
        "latex", help="Generate LaTeX table and bib file from a NNPDF runcard."
    )
    latex_parser.add_argument(
        "runcard",
        help="Path to the runcard .yaml file containing a ``dataset_inputs`` entry.",
        type=Path,
    )
    latex_parser.add_argument(
        "--group-by", action="store_true", help="Generate a separate latex table per group."
    )
    latex_parser.add_argument(
        "--output-bib",
        type=Path,
        help="Output file for the .bib file (if not given, prints to stdout)",
    )

    # Make the sorting option as args.sort_mode
    for option_parser in [latex_parser, list_parser]:
        sort_group = option_parser.add_mutually_exclusive_group()

        for option, help_text in [
            ("experiment", "Order by experiment."),
            ("process", "Order by process."),
            ("nnpdf31-process", "Order by NNPDF3.1 process."),
        ]:
            sort_group.add_argument(
                f"--{option}",
                dest="sort_mode",
                action="store_const",
                const=option.replace("-", "_"),
                help=help_text,
            )

    args = parser.parse_args()
    _setup_logging()

    # Now, regardless of the command, we need to build the whole index of dataset
    # NB: caching this would be nice but it doesn't seem to be too heavy anywhere
    # and the rules to invalidate the cache makes the rest of the code unnecessarily complicated.
    entries = index.build_index()

    if args.command == "list":
        _cmd_list(entries, args.filter, args.sort_mode, args.yaml)
    elif args.command == "latex":
        runcard_path = args.runcard

        if args.group_by:
            group_by = args.sort_mode
        else:
            group_by = None
        _cmd_latex(
            entries, runcard_path, args.sort_mode, group_by=group_by, output_bib=args.output_bib
        )
