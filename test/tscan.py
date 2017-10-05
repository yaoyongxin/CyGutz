#!/usr/bin/env python

import glob, os, subprocess

for path in glob.glob('test*'):
    os.chdir(path)
    cmd = ['python', path+'.py']
    subprocess.call(cmd)
