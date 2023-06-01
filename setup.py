"""Setup script for the MEA module."""

from setuptools import find_packages
from setuptools import setup

import mea

setup(name='mea',
      description='Package to use model-guided exploration of sequence space to diagnostic design guide RNAs for CRISPR-Cas13a',
      url="https://github.com/broadinstitute/mea-cas13",
      version=1.0,
      packages=find_packages(),
      package_data={
        "mea": ["utils/gan_data/*"],
      },
      scripts=[
        'design_guides.py'
      ]) 