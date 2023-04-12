"""
replica_selector.py

Tools for filtering replica sets based on criteria on the replicas.
"""
import logging
import shutil
import re


from reportengine.checks import make_argcheck, check
from reportengine.compat import yaml

from validphys.core import PDF
from validphys.renametools import rename_pdf
from validphys.utils import tempfile_cleaner

log = logging.getLogger(__name__)

def _fixup_new_replica(alphas_pdf: PDF, new_replica_file):
    """Helper function that takes in a
    :py:class:`validphys.core.PDF` object as well as
    the path to the central replica corresponding to the
    PDF and handles the writing of the alphas values
    to the header file.
    """
    alphas_mz = alphas_pdf.alphas_mz
    alphas_vals = alphas_pdf.alphas_vals
    with open(new_replica_file, 'rb') as in_stream:
        data = in_stream.read()
    with open(new_replica_file, 'wb') as out_stream:
        # Add the AlphaS_MZ and AlphaS_Vals keys
        out_stream.write(f"AlphaS_MZ: {alphas_mz}\nAlphaS_Vals: {alphas_vals}\n".encode())
        out_stream.write(data)

@make_argcheck
def _check_target_name(target_name):
    """Make sure this specifies a name and not some kid of path"""
    if target_name is None:
        return
    check(
        re.fullmatch(r'[\w]+', target_name),
        "`target_name` must contain alphnumeric characters and underscores only",
    )

@_check_target_name
def alpha_s_bundle_pdf(pdf, pdfs, output_path, target_name: (str, type(None)) = None):
    """Action that bundles PDFs for distributing to the LHAPDF
    format. The baseline pdf is declared as the ``pdf`` key
    and the PDFs from which the replica 0s are to be added is
    declared as the ``pdfs`` list.

    The bundled PDF set is stored inside the ``output`` directory.

    Parameters
    ----------
    pdf: :py:class:`validphys.core.PDF`
        The baseline PDF to which the new replicas will be added
    pdfs: list of :py:class:`validphys.core.PDF`
        The list of PDFs from which replica0 will be appended
    target_name: str or None
        Optional argument specifying the name of the output PDF.
        If ``None``, then the name of the original pdf is used
        but with ``_pdfas`` appended
    """
    base_pdf_path = pdf.infopath.parent
    nrep = len(pdf)

    target_name = target_name or pdf.name + '_pdfas'
    target_path = output_path / target_name

    alphas_paths = [i.infopath.parent for i in pdfs]
    alphas_replica0s = [path / f'{p}_0000.dat' for path, p in zip(alphas_paths, pdfs)]
    new_nrep = nrep + len(alphas_replica0s)
    alphas_values = [str(p.alphas_mz) for p in pdfs]

    if target_path.exists():
        log.warning(f"{target_path} already exists. Deleting contents.")
        shutil.rmtree(target_path)

    # We create a temporary directory to handle the manipulations inside.
    # We move the files to the new directory at the end.
    with tempfile_cleaner(
        root=output_path, exit_func=shutil.rmtree, exc=KeyboardInterrupt
    ) as tempdir:
        # Copy the base pdf into the temporary directory
        temp_pdf = shutil.copytree(base_pdf_path, tempdir / pdf.name)

        # Copy the alphas PDF replica0s into the new PDF
        for i, (alphas_pdf, rep) in enumerate(zip(pdfs, alphas_replica0s)):
            to = temp_pdf / f'{pdf.name}_{str(i + nrep).zfill(4)}.dat'
            shutil.copy(rep, to)
            _fixup_new_replica(alphas_pdf, to)

        #  Fixup the info file
        info_file = (temp_pdf / temp_pdf.name).with_suffix('.info')

        with open(info_file, 'r') as stream:
            yaml_obj = yaml.YAML()
            info_yaml = yaml_obj.load(stream)
        info_yaml['NumMembers'] = new_nrep
        info_yaml['ErrorType'] += '+as'
        extra_desc = '; '.join(
            f"mem={i} => alphas(MZ)={val}"
            for val, i in zip(alphas_values, range(nrep, new_nrep))
        )
        info_yaml['SetDesc'] += f"; {extra_desc}"
        with open(info_file, 'w') as stream:
            yaml_obj.dump(info_yaml, stream)

        # Rename the base pdf to the final name
        rename_pdf(temp_pdf, pdf.name, target_name)
        # This is the pdf path after the above renaming
        # i.e new_pdf.exists() == True
        new_pdf = temp_pdf.with_name(target_name)
        # Move the final pdf outside the temporary directory
        new_pdf = new_pdf.rename(target_path)
    log.info(f"alpha_s bundle written at {new_pdf}")
    return target_name

