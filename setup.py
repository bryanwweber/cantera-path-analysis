from setuptools import setup

setup(
    name='Cantera Path Analysis',
    version='0.0.1',
    packages=['cantera-path-analysis'],
    entry_points={
      'console_scripts': [
          'cantera-path-analysis = cantera-path-analysis.__main__:main'
      ]
    },
)
