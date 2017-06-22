# MRS_T2fit
Determine T2 relaxation times of MRS metabolites

The Software takes a Single Voxel Magnetic Resonance Spectroscopy 
dataset that contains multiple acquisitions at different echo times,
fits the Brain metabolites NAA, Cr, Cho and Glx and calculates their
T2 relaxation times.</br>
Knowledge of T2 is a pre-requisite for attenuation correction in
quantitative MR spectroscopy. Find a very useful introduction to the
subject at:  https://pgmm03.wordpress.com/2015/05/01/relax-dont-do-it/

This software is basically a script calling the MRS processing tool</br>
TARQUIN (http://tarquin.sourceforge.net) 

### Requirements
&emsp;<a href="http://www.python.org">- Python</a></br>
&emsp;<a href="http://www.numpy.org/">- Numpy</a></br>
&emsp;<a href="http://www.scipy.org/">- SciPy</a></br>
The following additional Python libraries are strongly recommended:</br>
(the software works without, but with some loss of functionality)</br>
&emsp;<a href="http://pydicom.readthedocs.io">- PyDicom</a> (to read spectro files in DICOM format)</br> 
&emsp;<a href="http://wiki.python.org/moin/TkInter">- TkInter</a> (to interactively choose input files)</br> 
&emsp;<a href="http://pypi.python.org/pypi/pywin32">- PyWin32</a> (to gracefully handle user interrupts on Windows)</br>
The program also requires the following external files from <a href="http://tarquin.sourceforge.net">TARQUIN</a></br>
&emsp;- Windows: tarquin.exe, cvm_ia32.dll, libfftw-3.2.2.dll, vcomp100.dll</br>
&emsp;- MacOS:   tarquin, libz.1.dylib, libpng15.15.dylib</br>
&emsp;- Linux:   tarquin</br>
For convenience these files are included in the "Supportfiles_" archives for the respective platforms.</br>
Extract and copy to the folder were the main python script resides.</br>

### Runs under:
- Windows (standalone binary included under the release tab)
- Linux
- MacOS

### Usage:
    MRS_T2fit.py --spec=<spectrofile>
    MRS_T2fit.py --help

### MR data:
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) <b> Currently supports Philips formats only </b> ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

- SPAR/SDAT format, or
- Philips DICOM format(spectra are in files called XX*)

The dataset MUST contain acquisitions of at least 7 different TE's acquired with a SV PRESS sequence.

For convenience a suggestion for acquisition parameters is included in "Acquisition Parameters.txt"

##
### License:
<a href="http://www.gnu.org/licenses">GPLv2</a> (see also inside `MRS_T2fit.py`)

### Author:
Bernd Foerster, bfoerster at gmail dot com
