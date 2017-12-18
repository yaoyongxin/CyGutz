#!/usr/bin/env python

import glob, os, subprocess

cwd = os.getcwd()
for path in glob.glob('test*'):
    os.chdir(cwd+'/'+path)
    if os.path.isfile(path+".py"):
        cmd = ["python", path+".py"]
        subprocess.call(cmd)
