from setuptools import setup 

with open("README.md", "r") as f:
    long_desc = f.read()

setup(
        name='pyarcknot', 
        version='0.0.2',
        description='A package for studying the arc diagrams of mathematical knots',
        long_description=long_desc,
        long_description_content_type="text/markdown",
        author='inigorrix', 
        author_email='inigorrix@gmail.com',
        url='https://github.com/inigorrix/pyarcknot',
        packages=['pyarcknot'],
        license="MIT",
        install_requires=[ 
	        'numpy',
	        'matplotlib',
	        'sympy',
	        'numba'],
        keywords=[
                'knot',
                'knot diagram',
                'arc diagram', 
                'jones polynomial',
                'turaev surface'],
        classifiers=[
                'License :: OSI Approved :: MIT License',
                'Operating System :: OS Independent',
                'Programming Language :: Python :: 3']
) 
