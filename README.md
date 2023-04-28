SPARC4 Quality Assurance Scripts
================================

This repository contains scripts to verify the images from the [SPARC4 - Simultaneous Polarimeter and Rapid Camera in Four bands](http://www.das.inpe.br/sparc4/), installed at [Observatório Pico dos Dias - Brazil](https://www.gov.br/lna/pt-br/composicao-1/coast/obs/opd) (OPD/LNA/MCTI).

All scripts are meant to be used internally by the SPARC4 Developed team.

Install
-------

This repository do not contain a python package by itself. It just contain scripts. So, you don't have to install it. Just download the scripts, install the dependencies of each one, and execute the scripts using python:

```
python <script.py> <options...>
```

Scripts
-------

- **sparc4_raw_header_check.py**: check if all the header keywords comply with the standard stabilished during instrument and data design. Any keyword that doesn't match the standard is printed, with the error and the standard. *This script depends on astropy package.*

```bash
# to simply print in terminal
python sparc4_raw_header_check.py file1.fits file2.fits ...
# to save output to txt
python sparc4_raw_header_check.py -o report.txt -f txt file1.fits file2.fits ...
# to save output to csv
python sparc4_raw_header_check.py -o report.csv -f csv file1.fits file2.fits ...
# to print changelog
python sparc4_raw_header_check.py --changelog
# to print version
python sparc4_raw_header_check.py --version
# to output the standard to a csv table
python sparc4_raw_header_check.py --print-standard standard.csv
```
