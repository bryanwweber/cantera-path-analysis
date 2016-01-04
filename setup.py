from setuptools import setup

setup(
    name='pychempath',
    version='0.0.1',
    packages=['pychempath'],
    entry_points={
      'console_scripts': [
          'pychempath = pychempath.__main__:main',
      ],
    },
)
