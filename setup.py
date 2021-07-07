# -*- coding: utf-8 -*-
from setuptools import setup

packages = \
['ibs']

package_data = \
{'': ['*']}

install_requires = \
['Jinja2>=3.0.1,<4.0.0', 'click>=7.0.0,<8.0.0', 'matplotlib>=3.4.2,<4.0.0']

entry_points = \
{'console_scripts': ['runode = ibs:cli_runode.main']}

setup_kwargs = {
    'name': 'ibs',
    'version': '0.0.0',
    'description': '<Enter a one-sentence description of this project here.>',
    'long_description': '===========================\nIBS - Intra Beam Scattering\n===========================\n\n\nC++ Library for IBS with Python wrapper.\n\n\n* Free software: MIT license\n* Documentation: https://ibs.readthedocs.io.\n\n\nFeatures\n--------\n\n* Twiss module\n',
    'author': 'Tom Mertens',
    'author_email': 'your.email@whatev.er',
    'maintainer': None,
    'maintainer_email': None,
    'url': 'https://github.com/tomerten/ibs',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'entry_points': entry_points,
    'python_requires': '>=3.7,<4.0',
}
from build import *
build(setup_kwargs)

setup(**setup_kwargs)
