#We don't want to import stuff here that could slow down the import times
def __give_git():
    from subprocess import run
    return run('git describe --tags --dirty'.split(), capture_output=True).stdout.decode().strip()
__build__ = __give_git()
