# Copyright 2018 The Cirq Developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import io
from setuptools import setup

# This reads the __version__ variable from version.py
__version__ = ''
exec(open('chp_sim/version.py').read())

description = ("Reference implementation of Aaronson et al's CHP simulator "
               "for efficiently simulating quantum stabilizer circuits.")

# README file as long_description.
long_description = io.open('README.md', encoding='utf-8').read()

# Read in requirements
requirements = open('requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setup(name='chp_sim',
      version=__version__,
      url='https://github.com/Strilanc/python-chp-stabilizer-simulator',
      author='Craig Gidney',
      author_email='craig.gidney@gmail.com',
      python_requires='>=3.6.0',
      install_requires=requirements,
      license='Apache 2',
      description='',
      long_description=long_description,
      packages=['chp_sim'])
