from setuptools import setup, find_packages


scm_version_template = """# Generated by setuptools_scm
__all__ = ["__version__"]

__version__ = "{version}"
"""

setup(
    name="syseng_throughputs",
    use_scm_version={
        "write_to": "syseng_throughputs/version.py",
        "write_to_template": scm_version_template,
    },
    scripts=[
      "bin/calcM5"
    ],
    packages=find_packages(),
)
