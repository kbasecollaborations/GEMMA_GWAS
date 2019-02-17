#!/usr/bin/python

import os
import subprocess

files = os.listdir('.')

for file in files:
    if file.endswith('.vcf'):
        cmd = ['vcf_validator', '-i', file]
        subprocess.call(cmd)
