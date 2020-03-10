"""
    A collection of utility functions to handle logistics of
    LHAPDFs and fits. For use by vp-scripts.
"""
import os

def rename_pdf(pdf_folder, initial_fit_name, final_name):
    for item in os.listdir(pdf_folder):
        p = pdf_folder/item
        if p.is_symlink():
            replica = p.resolve().parent.name
            pointer = f'../../nnfit/{replica}/{final_name}.dat'
            p.unlink()
            p.symlink_to(pointer)
        newname = p.name.replace(initial_fit_name, final_name)
        p.rename(p.with_name(newname))
    pdf_folder.rename(pdf_folder.with_name(final_name))


def rename_nnfit(nnfit_path, initial_fit_name, final_name):
    info_file = nnfit_path/f'{initial_fit_name}.info'
    info_file.rename(info_file.with_name(f'{final_name}.info'))
    #Some older fits have the PDF here
    pdf_folder = nnfit_path / initial_fit_name
    if pdf_folder.is_dir():
        rename_pdf(pdf_folder, initial_fit_name, final_name)
    #Change replica names
    for item in nnfit_path.glob('replica*'):
        if item.is_dir():
            files = item.glob(initial_fit_name + '*')
            for i in files:
                newname = i.name.replace(initial_fit_name, final_name)
                i.rename(item/newname)


def rename_postfit(postfit_path, initial_fit_name, final_name):
    pdf_folder = postfit_path / initial_fit_name
    rename_pdf(pdf_folder, initial_fit_name, final_name)
    os.system(f'sed -i -e "s/{initial_fit_name}/{final_name}/g" {postfit_path/"postfit.log"}')

def change_name(initial_path, final_name):
    """Function that takes initial fit name and final fit name
    and performs the renaming"""
    initial_fit_name = initial_path.name
    nnfit = initial_path/'nnfit'
    if nnfit.exists():
        rename_nnfit(nnfit, initial_fit_name, final_name)
    postfit = initial_path/'postfit'
    if postfit.exists():
        rename_postfit(postfit, initial_fit_name, final_name)
    newpath = initial_path.with_name(final_name)
    initial_path.rename(newpath)
    return newpath

