import os

DIR = os.path.abspath(os.path.dirname(__file__))


def get_include() -> str:
    installed_path = os.path.join(DIR, "include")
    source_path = os.path.join(os.path.dirname(DIR), "include")
    return installed_path if os.path.exists(installed_path) else source_path
