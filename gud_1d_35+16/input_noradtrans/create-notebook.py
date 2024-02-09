#!/usr/bin/env python3
"""create-notebook.py
   Creates a minimal jupyter notebook (.ipynb)
   Usage: ./create-notebook.py <notebook>

   ref: https://stackoverflow.com/questions/48946639/how-to-create-and-open-a-jupyter-notebook-ipynb-file-directly-from-terminal
   answer by stacksonstacks
"""
import sys
from notebook import transutils as _
from notebook.services.contents.filemanager import FileContentsManager as FCM

try:
    notebook_fname = sys.argv[1].strip('.ipynb')
except IndexError:
    print("Usage: create-notebook <notebook>")
    exit()

notebook_fname += '.ipynb'  # ensure .ipynb suffix is added
FCM().new(path=notebook_fname)
