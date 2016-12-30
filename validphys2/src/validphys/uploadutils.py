#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 10:38:21 2016

@author: zah
"""
import subprocess
import logging
import shutil
import uuid
import base64

from reportengine.colors import t

log = logging.getLogger(__name__)

#TODO: This should go to nnprofile
upload_host = 'nnpdf@pcteserver.mi.infn.it'

target_dir = "WEB/validphys-reports/"

root_url = 'http://pcteserver.mi.infn.it/~nnpdf/validphys-reports/'

class BadSSH(Exception): pass

def check_auth():
    ssh_command_line = ('ssh', '-o', 'PreferredAuthentications=publickey',
                        '-q', upload_host, 'exit')

    str_line = ' '.join(repr(ele) for ele in ssh_command_line)

    log.info("Checking SSH connection to %s.", upload_host)

    try:
        subprocess.run(ssh_command_line, check=True)
    except subprocess.CalledProcessError as e:
        raise BadSSH(("Could not validate the SSH key. "
        "The command\n%s\nreturned a non zero exit status. "
        "Please make sure thet your public SSH key is on the server.")
        % str_line) from e
    except OSError as e:
        raise BadSSH("Could not run the command\n%s\n: %s" % (str_line, e)) from e

    log.info("Connection seems OK.")

def check_rsync():
    if not shutil.which('rsync'):
        raise BadSSH("Could not find the rsync command. "
        "Please make sure it is installed.")

def check_upload():
    check_rsync()
    check_auth()


def upload_output(output_path):
    randname = base64.urlsafe_b64encode(uuid.uuid4().bytes).decode()
    newdir = target_dir + randname

    rsync_command = ('rsync', '-az',
                     str(output_path)+'/',
                     ':'.join((upload_host, newdir)))

    log.info("Uploading output (%s) to %s" % (output_path, upload_host))
    try:
        subprocess.run(rsync_command, check=True)
    except subprocess.CalledProcessError as e:
        msg = "Failed to upload output: %s" % e
        raise BadSSH(msg) from e
    else:
        url = root_url + randname
        log.info("Upload completed. The result is available at:\n%s" % t.bold_blue(url))



