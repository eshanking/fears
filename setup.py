

from setuptools import setup, find_packages

setup(
    name='fears', 
    version='0.1.0', 
    author = 'Eshan King', 
    author_email = '',
    # packages=['autorate', "autorate.test", "autorate.test.data"], 
    # pacakges=['fears','fears.data','fears.utils'],
    packages=find_packages(include=['fears','fears.population','fears.data','fears.utils','fears.utils.*']),
    # packages=find_packages(where="fears"),
    install_requires = [
      "pandas",
      "pytest",
      "scipy",
      "matplotlib",
      "numpy",
      "importlib_resources", 
      "lifelines", 
      "seaborn", 
      "networkx",
      "cycler",
      "matplotlib-label-lines"
    ],
    include_package_data=True,
    package_data={'': ['data/*.csv']}
)