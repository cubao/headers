from setuptools import setup
import os
import shutil


from tempfile import TemporaryDirectory

with TemporaryDirectory() as temp_dir:
    base_dir = os.path.abspath(os.path.dirname(__file__))

    for name in ["include"]:
        shutil.copytree(
            os.path.join(base_dir, name),
            os.path.join(temp_dir, name),
            dirs_exist_ok=True,
        )

    setup(
        name="cubao_headers",
        version="0.0.1",
        author="district10",
        author_email="dvorak4tzx@gmail.com",
        description="Universal headers used by cubao team.",
        url="https://github.com/cubao/headers",
        license="MIT",
        packages=["cubao_headers"],
        zip_safe=False,
        package_dir={"cubao_headers": temp_dir},
        package_data={
            "cubao_headers": [
                "include/**/*",
            ]
        },
    )
