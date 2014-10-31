from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()
    
setup(name='nsga2lib',
      version='2.0',
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
      include_package_data=True
      )
