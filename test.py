import os
import glob
from pathlib import Path

os.chdir('include')
dirs = sorted(glob.glob('**/*/', recursive=True))
for dir in dirs:
    line = 'f"{name}.include.' + dir[:-1].replace('/', '.') + '": pattern,'
    print(line)
    path = f'{dir}/__init__.py'
    if os.path.isfile(path):
        continue
    Path(path).touch()
