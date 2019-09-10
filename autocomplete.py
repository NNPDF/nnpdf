#!/usr/bin/env python3
import sys

from importlib import import_module
from os import path

from reportengine.utils import get_providers
from validphys.app import providers


def candidates():
    if not path.exists("./providers.txt"):
        dictionary = [get_providers(import_module(m)) for m in providers]
        functions = [item for sublist in [list(i.keys()) for i in dictionary]
                    for item in sublist]
        with open("providers.txt", "w") as stream:
            for function in functions:
                stream.write(f"{function}\n")
    else:
        with open("providers.txt", "r") as stream:
            functions = stream.read().split("\n")
    return functions


def completion_hook(cmd, curr_word, prev_word):
    potential_matches = candidates()
    matches = [k for k in potential_matches if k.startswith(curr_word)]
    return matches


def main():
    if sys.argv[-1] not in ["-h", "--help"]:
        return

    results = completion_hook(*sys.argv[1:])
    if len(results):
        print("\n".join(results))


if __name__ == "__main__":
    main()
