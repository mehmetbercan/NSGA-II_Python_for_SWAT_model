from setuptools import setup
import os

def readme():
    with open('README.rst') as f:
        return f.read()
    
setup(name='nsga2lib',
      version='2.2',
      description='Libraries for performing nsga2 calibration.',
      long_description=readme(),
      classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: Beta"
        "Intended Audience :: Science/Research",
        "Natural Language :: English",     
      ],
      url='https://github.com/mehmetbercan/NSGA-II_Python_for_SWAT_model',
      author='Mehmet B. Ercan',
      author_email='ercanm@email.sc.edu',
      license='MIT',
      packages=['nsga2lib'],
      install_requires=['nsga2lib','numpy'],
	  data_files   = [(os.path.join('nsga2lib','ScriptsForSWATtxt'),[os.path.join('nsga2lib','ScriptsForSWATtxt','Extract_rch.py'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','SWAT_ParameterEdit.py'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','Makefile'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','swat2012_627'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','nsga2_mid.sh'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','nsga2_mid.cmd'),
								  os.path.join('nsga2lib','ScriptsForSWATtxt','swat.exe'),])],
      include_package_data=True
      )
