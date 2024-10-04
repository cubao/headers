import os
import glob
from pathlib import Path

head = '''#!/usr/bin/env python3
from setuptools import setup, find_packages

name = "cubao_headers"
version = "0.0.7"

pattern = ["*"]

setup(
    name=name,
    version=version,'''
print(head)

os.chdir('cubao_headers/include')
path = f'__init__.py'
Path(path).touch()

lines = '''import os

DIR = os.path.abspath(os.path.dirname(__file__))


def get_include() -> str:
    installed_path = os.path.join(DIR, "include")
    source_path = os.path.join(os.path.dirname(DIR), "include")
    return installed_path if os.path.exists(installed_path) else source_path
'''

with open('../__init__.py', 'w') as f:
    f.write(lines)


packages = ['f"{name}"', 'f"{name}.include"']
prefix = '        '

dirs = sorted(glob.glob('**/*/', recursive=True))
lines = []
for dir in dirs:
    package = 'f"{name}.include.' + dir[:-1].replace('/', '.') + '"'
    packages.append(package)
    line = f'{prefix}{package}: pattern,'
    lines.append(line)
    path = f'{dir}__init__.py'
    Path(path).touch()

print('    packages=[' + ',\n        '.join(packages) + '],')
print('    package_data={')
line = 'f"{name}.include": pattern,'
print(f'{prefix}{line}')
for line in lines:
    print(line)
tail = '''    },
)'''
print(tail)
