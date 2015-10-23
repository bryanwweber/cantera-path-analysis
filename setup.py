from setuptools import setup

setup(
    name='CanPathAn',
    version='0.0.1',
    packages=['canpathan'],
    entry_points={
      'console_scripts': [
          'cantera-path-analysis = canpathan.__main__:main',
      ],
    },
)
