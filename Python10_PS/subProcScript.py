#!/usr/bin/env	python3
import sys
import subprocess

output = subprocess.check_output('ls -l | grep change' , shell = True)
print(output.decode('utf-8')) 
