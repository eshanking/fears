from setuptools import setup

setup(
    name='fears', 
    version='0.9.0', 
    author = 'Eshan King', 
    author_email = 'eshan.king@case.edu',
    packages=['fears', 'fears.classes', 'fears.utils'], 
    install_requires = [
      "pathlib", 
      "pandas",
      "numpy",
      "pickle",
      "lifelines", 
      "networkx", 
      "matplotlib", 
      "cycler",
      "seaborn",
      "scipy", 
      "labellines" 
    ],
    python_requires='>= 3.6.2, <=3.9.7'
)