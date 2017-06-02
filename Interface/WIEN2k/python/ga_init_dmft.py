#!/usr/bin/env python

import optparse
import os, shutil
import subprocess
import ga_indmffile
from ga_utils import DmftEnvironment, W2kEnvironment


def print_prompt(text):
    print "----->", text


dmfenv = DmftEnvironment()
w2kenv = W2kEnvironment()

indmf = ga_indmffile.Indmf(w2kenv.case)

# First, check if case.indmf exists
create_indmf = True
if indmf.file_exists():
    prompt = "File `%s` already exists.  Do you want to create new inputs? (y/n): " % indmf.filename()
    userin = raw_input(prompt).lower().strip()
    if userin.lower() == 'n':
        create_indmf = False

# create new case.indmf (if desired)
if create_indmf:
    indmf.user_input()
    indmf.write()
    fname = indmf.filename()
    fnametemp = fname + '_st'

    finished_editing = False
    while not finished_editing:
        shutil.copy(fname, fnametemp)  # copy to temp file for user to edit
        print_prompt("Edit %s file." % (fnametemp,))
        subprocess.call(w2kenv.EDITOR.split() + [fnametemp])

        # read user-edited file in just to check for any syntax errors
        indmf_st = ga_indmffile.Indmf(w2kenv.case)
        indmf_st.extn = 'indmf_st'
        try:
            indmf_st.read()
            finished_editing = True
            shutil.copy(fnametemp, fname)  # move it back            
        except:
            print_prompt("Your edits have introduced syntax errors in %s." % (fnametemp,))
            prompt = "Do you want to edit the file again? (y/n): "
            userin = raw_input(prompt).lower().strip()
            if userin.lower() == 'n':
                finished_editing = True

# ask user if this is a SO run
so_run = False
prompt = "Is this a spin-orbit run? (y/n): "
userin = raw_input(prompt).lower().strip()
if userin.lower() == 'y':
    so_run = True

# run sigen to calculate siginds, cftrans and write output to case.indmfl file
sigen = os.path.join( dmfenv.ROOT, 'ga_sigen.py' )
# sigen = './sigen.py'

print sigen

args = ["--so"] if so_run else []
if subprocess.call([sigen] + args):
    print "Error executing sigen.py."
    exit(1)

# copy case.indmfl --> case.indmfl_st to allow user to make changes
findmfl = w2kenv.case + '.indmfl'
findmfltemp = findmfl + '_st'

finished_editing = False
while not finished_editing:
    shutil.copy(findmfl, findmfltemp)
    print_prompt("Edit %s file." % (findmfltemp,))
    subprocess.call(w2kenv.EDITOR.split() + [findmfltemp])

    # read user-edited file in just to check for any syntax errors
    inl = ga_indmffile.Indmfl(w2kenv.case)
    inl.extn = 'indmfl_st'
    try:
        inl.read()
        finished_editing = True
        shutil.copy(findmfltemp, findmfl)  # move it back
    except:
        print_prompt("Your edits have introduced syntax errors in %s." % (findmfltemp,))
        prompt = "Do you want to edit the file again? (y/n): "
        userin = raw_input(prompt).lower().strip()
        if userin.lower() == 'n':
            finished_editing = True

