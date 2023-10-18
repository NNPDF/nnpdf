# This file gets overwritten on deployemt, don't modify
def __give_git():
    # Get the path to this file
    from pathlib import Path

    file_dir = Path(__file__).parent
    from subprocess import CalledProcessError, run

    try:
        result = run(
            ["git", "describe", "--tags", "--long", "--dirty=-dev"],
            capture_output=True,
            text=True,
            check=True,
            cwd=file_dir,
        ).stdout.strip()
        tag = run(
            ["git", "describe", "--abbrev=0", "--tags"],
            capture_output=True,
            text=True,
            check=True,
            cwd=file_dir,
        ).stdout.strip()
        githash = run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=file_dir,
        ).stdout.strip()
        version = result.replace(f"-g{githash}", f"+g{githash}").replace(f"{tag}-", f"{tag}.")
    except CalledProcessError:
        # In principle this function should not exist on an installed version
        # but who knows. Also maybe git doesn't work on the machine or whatever
        version = "unknown"
    return version


build_version = __give_git()
