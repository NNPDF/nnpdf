# This file gets overwritten on deployemt, don't modify
def __give_git():
    # Get the path to this file
    from pathlib import Path
    file_dir = Path(__file__).parent
    from subprocess import run, CalledProcessError
    try:
        result = run('git describe --tags --dirty'.split(), capture_output=True, text=True, check=True, cwd=file_dir).stdout.strip()
    except CalledProcessError:
        # In principle this function should not exist on an installed version
        # but who knows. Also maybe git doesn't work on the machine or whatever
        result = 'unknown'
    return result
build_version = __give_git()
