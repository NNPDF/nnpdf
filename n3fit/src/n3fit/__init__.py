def __give_git():
    from subprocess import run
    return run('git describe --tags --dirty'.split(), capture_output=True).stdout.decode().strip()
__build__ = __give_git()
