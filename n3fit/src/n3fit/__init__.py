def __give_version():
    # Get the path to this file
    from pathlib import Path
    file_dir = Path(__file__).parent
    from subprocess import run, CalledProcessError
    try:
        result = run('git describe --tags --dirty'.split(), capture_output=True, text=True, check=True, cwd=file_dir).stdout.strip()
    except CalledProcessError:
        from pkg_resources import get_distribution
        # If we cannot get the name directly form git for whatever the reason, ask the package
        result = get_distribution(__name__).version
    return result
__version__ = __give_version()
