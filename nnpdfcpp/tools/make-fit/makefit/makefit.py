# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:54:18 2016

@author: Zahari Kassabov
"""
import argparse
import pathlib
import shutil
import subprocess
import datetime
import sys

import yaml
import blessings

CONFIG_VERSION = '1.0'

nnconfig = pathlib.Path('~/.nnpdfrc').expanduser()
config_keys = ('user', 'initials', 'editor', 'build_path')

def yes_or_no(prompt, default = None):
    yes = {'y', 'yes'}
    no = {'n', 'no'}
    if default is None:
        appendix = "[y/n]"
    elif default:
        appendix = "[" + t.bold('Y') + "/n]:\n"
    else:
        appendix = "[y/"+t.bold('N')+"]:\n"
    while True:
        try:
            val = input(prompt + " " + appendix)
        except (EOFError, KeyboardInterrupt):
            sys.exit(1)
        if val == '' and default is not None:
            return default
        if val.lower() in yes:
            return True
        elif val.lower() in no:
            return False
        print("Please answer 'yes' or 'no'")

#NOT USED ANYMORE
def get_git_rev(path):
    path = str(path)
    call = subprocess.run(['git', 'rev-parse', 'HEAD'],
                   stdout=subprocess.PIPE, universal_newlines=True,
                   cwd=path, check=True)
    return call.stdout[:7]


t = blessings.Terminal()

def ask_for_user_and_initials():
    while True:
        user = input("Please enter your " + t.bold("name and surname:\n"))
        bits = user.split()
        if len(bits) > 1:
            initials = ''.join(x[0].lower() for x in bits)
            print("Your initials are '%s'" % t.blue(initials))
            break
        else:
            print("Please enter at least " +  t.bold_red("two words"))

    return user, initials

def ask_for_editor():
    while True:
        editor = input("Please enter your favourite " + t.bold('editor') +
                       " (for yaml files) [default: vi]:\n")
        if editor == '':
            editor = 'vi'
        editor_path = shutil.which(editor)
        if editor_path:
            print("Found editor: %s" % t.blue(editor_path))
            break
        else:
            if yes_or_no("Editor '%s' not found in path. Are you sure" % editor, False):
                break
    return editor

#TODO: Trash all this and just use the conda build str someday.
def ask_for_path():
    here = pathlib.Path('.').resolve()
    if here.stem == 'nnpdfcpp':
        default = here / 'nnpdfbuild'
        print(default)
    else:
        try:
            default = next(x for x in here.parents if x.stem == 'nnpdfbuild')
        except StopIteration:
            default = None

    prompt =  "Please enter the path to the 'nnpdfbuild' directory"
    if default:
        prompt += " [default: %s]" % default
        prompt += ':'
    while True:
        val = input(prompt)
        if val == '' and default:
            build_path = default
            print("build_path = default")
        else:
            build_path = val
        try:
            print(build_path)
            p = pathlib.Path(build_path).resolve()
            print(t.yellow(str(p)))
        except OSError as e:
            print(t.red_bold("Bad path: %s" %e))
            continue
        try:
            rev = get_git_rev(p)
        except Exception as e:
            print(t.red_bold("Could not get git rev fro path %s" %e ))
            continue
        print("Current git rev is %s" % t.blue(rev))
        print(p)
        return str(p)

def ask_for_config():
    user, initials = ask_for_user_and_initials()
    editor = ask_for_editor()
    build_path = ask_for_path()

    with nnconfig.open('w') as f:
        yaml.dump({'user': user,'initials': initials,
                   'editor': editor,
                   'build_path': build_path},
                   f)
    print("Settings saved in %s" % t.blue(str(nnconfig)))

def get_config():
    try:
        with nnconfig.open() as f:
            config = yaml.load(f)
    except FileNotFoundError:
        ask_for_config()
        return get_config()

    if not isinstance(config, dict):
        print(t.red_bold("Bad config. Regenerating"))
        ask_for_config()
        return get_config()

    ret = {}

    for key in config_keys:
        if not key in config:
            print(t.red_bold("Key '%s' missing in config." % key) +
                 " Generating config again.")
            ask_for_config()
            get_config()
            return get_config()

        ret[key] = config[key]
    return ret

def make_name(config):
    config_path = pathlib.Path(config['build_path']) / 'config'
    initials = config['initials']
    now = datetime.datetime.now()
    now = format(now, "%g%m%d")

    i = 1
    while True:
        path = config_path / ('-'.join((now, '%03d'%i ,initials)) + '.yml')
        if not path.exists():
            break
        i+=1
    return str(path)


def main():
    config = get_config()

    parser = argparse.ArgumentParser(description="Generates a new fit config "
                                                 "from a reference")

    parser.add_argument("reference")

    args = parser.parse_args()

    try:
        with open(args.reference) as f:
            reference_text = f.read()
        reference_content = yaml.load(reference_text)
    except Exception as e:
            print(t.red_bold("Could not load refence file: %s" %e))
            sys.exit()

    try:
        name = make_name(config)
    except Exception as e:
        print(t.red_bold("Could not make name: %s" % e))
        raise e
        sys.exit(1)
    shutil.copy(args.reference,name)
    while True:
        subprocess.run([config['editor'], name])
        subprocess.run(['yamldiff', args.reference, name,
                        '--set-keys', 'experiments:experiment', 'datasets:dataset'])
        if yes_or_no("Do you accept the changes?"):
            break
    print(t.bold("Written new config file %s" % name))

if __name__ == '__main__':
    main()
