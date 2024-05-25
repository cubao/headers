import os
import glob
from pathlib import Path

head = '''#!/usr/bin/env python3
from setuptools import setup, find_packages

name = "cubao_headers"
version = "0.0.3"

pattern = ["*"]

setup(
    name=name,
    version=version,
    packages=find_packages(),
    package_data={'''
print(head)

os.chdir('src/cubao_headers/include')
path = f'__init__.py'
Path(path).touch()
path = f'../__init__.py'
Path(path).touch()

prefix = '        '
line = 'f"{name}.include": pattern,'
print(f'{prefix}{line}')
dirs = sorted(glob.glob('**/*/', recursive=True))
for dir in dirs:
    line = 'f"{name}.include.' + dir[:-1].replace('/', '.') + '": pattern,'
    print(f'{prefix}{line}')
    path = f'{dir}__init__.py'
    Path(path).touch()
tail = '''    },
)'''
print(tail)
