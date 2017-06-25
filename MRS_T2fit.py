#!/usr/bin/python
#
# MRS_T2fit - Determine T2 relaxation times of MRS metabolites
#     
# author: Bernd Foerster, bfoerster at gmail dot com
# 
# ----- VERSION HISTORY -----
#
# Version 0.2 - 22, June 2017
#   - public release on GitHub
#    
# ----- LICENSE -----                 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License (GPL) as published 
# by the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version. For more detail see the 
# GNU General Public License at <http://www.gnu.org/licenses/>.
# also the following additional license terms apply:
#
# ----- REQUIREMENTS ----- 
#
#   This program was developed under Python Version 2.7.6 for Windows 32bit
#   A windows standalone executable can be "compiled" with PyInstaller
#
#   The following Python libraries are required:
#     - NumPy (http://www.numpy.org/)
#     - SciPy (http://www.scipy.org/)
#   The following additional Python libraries are strongly recommended:
#   (the software works without, but with some loss of functionality)
#     To read spectro files in DICOM format 
#     - pydicom (http://pydicom.readthedocs.io)
#     To interactively choose input files 
#     - tkinter (http://wiki.python.org/moin/TkInter)
#     To gracefully handle user interrupts on Windows
#     - pywin32 (http://pypi.python.org/pypi/pywin32) 
#     
#   The program also requires the following external files
#   from TARQUIN (http://tarquin.sourceforge.net/):
#     - Windows: tarquin.exe, cvm_ia32.dll, libfftw-3.2.2.dll, vcomp100.dll
#     - MacOS:   tarquin, libz.1.dylib, libpng15.15.dylib
#     - Linux:   tarquin
#
# ----- TO DO LIST ----- 
#
#   - accept files from other vendors (currently Philips only)
#


Program_version = "v0.2" # program version


import sys
import math
import os
import signal
import random
import shutil
import subprocess
import time
import datetime
from getopt import getopt
from getopt import GetoptError
from distutils.version import LooseVersion

import csv
import numpy
from scipy.optimize import curve_fit
try: 
    from scipy.sparse.csgraph import _validation    # needed for pyinstaller
    from scipy.special import _ufuncs_cxx           # needed for pyinstaller
except: pass

TK_installed=True
try: from tkFileDialog import askopenfilename # Python 2
except: 
  try: from tkinter.filedialog import askopenfilename; # Python3
  except: TK_installed=False
try: import Tkinter as tk; # Python2
except: 
  try: import tkinter as tk; # Python3
  except: TK_installed=False

FNULL = open(os.devnull, 'w')
old_target, sys.stderr = sys.stderr, FNULL # replace sys.stdout 
pydicom_installed=True
try: import dicom
except: pydicom_installed=False
sys.stderr = old_target # re-enable

pywin32_installed=True
try: import win32console, win32gui, win32con
except: pywin32_installed=True


TE1 = 17 #[ms] - time between the 90' excitation pulse and the first PRESS echo
dTE = 10 #[ms] - TE step (only used when reading SPAR, for DICOM the actual TE values are used)
metabolites = numpy.array(['TNAA', 'TCho', 'TCr', 'Glx']) # only these will be evaluated
min_nTE=7 # minimum number of TE's required
if sys.platform=="win32": slash='\\'
else: slash='/'

def exit (code):
    # cleanup 
    try: shutil.rmtree(tempdir)
    except: pass # silent
    if pywin32_installed:
        try: # reenable console windows close button (useful if called command line or batch file)
            hwnd = win32console.GetConsoleWindow()
            hMenu = win32gui.GetSystemMenu(hwnd, 1)
            win32gui.DeleteMenu(hMenu, win32con.SC_CLOSE, win32con.MF_BYCOMMAND)
        except: pass #silent
    sys.exit(code)
def signal_handler(signal, frame):
    lprint ('User abort')
    exit(1)
def logwrite(message): 
    sys.stderr.write(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    sys.stderr.write(' ('+ID+') - '+message+'\n')
    sys.stderr.flush()   
def lprint (message):
    print (message)
    logwrite(message)
def checkfile(file): # generic check if file exists
    if not os.path.isfile(file): 
        lprint ('ERROR:  File "'+file+'" not found '); exit(1)    
def expdecay(x,A,T2): # T2 of long component fixed to CSF_T2
    return A*numpy.exp(-x/T2)
def isDICOM (file):
    try: f = open(file, "rb")
    except: lprint ('ERROR: opening file'+file); exit(1)
    try:
        test = f.read(128) # through the first 128 bytes away
        test = f.read(4) # this should be "DICM"
        f.close()
    except: return False # on error probably not a DICOM file
    if test == "DICM": return True 
    else: return False    
def _get_from_SPAR (input, varstring):
    if varstring[len(varstring)-1] != ' ': varstring += ' ' # requires final space
    value = [text.split(':')[1] for text in input if text.split(':')[0]==varstring]
    if len(value)>0: value = value[0]
    else: lprint ('ERROR: unable to read parameter "'+varstring+'" in SPAR'); sysexit(1)
    return value        
def delete (file):
    try: os.remove(file)
    except: pass #silent
def run (command, parameters):
    string = '"'+command+'" '+parameters
    if debug: logwrite (string)
    process = subprocess.Popen(string, env=my_env,
                  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = process.communicate()  
    if debug: logwrite (stdout)
    if debug: logwrite (stderr)    
    if process.returncode != 0: 
        lprint ('ERROR:  returned from "'+os.path.basename(command)+
                '", for details inspect logfile in debug mode')
        exit(1)
    return stdout    
def _vax_to_ieee_single_float(data): # borrowed from the python VeSPA project
    #Converts a float in Vax format to IEEE format.
    #data should be a single string of chars that have been read in from 
    #a binary file. These will be processed 4 at a time into float values.
    #Thus the total number of byte/chars in the string should be divisible
    #by 4.
    #Based on VAX data organization in a byte file, we need to do a bunch of 
    #bitwise operations to separate out the numbers that correspond to the
    #sign, the exponent and the fraction portions of this floating point
    #number
    #role :      S        EEEEEEEE      FFFFFFF      FFFFFFFF      FFFFFFFF
    #bits :      1        2      9      10                               32
    #bytes :     byte2           byte1               byte4         byte3    
    f = []; nfloat = int(len(data) / 4)
    for i in range(nfloat):
        byte2 = data[0 + i*4]; byte1 = data[1 + i*4]
        byte4 = data[2 + i*4]; byte3 = data[3 + i*4]
        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111
        sign  =  (ord(byte1) & 0x80) >> 7
        expon = ((ord(byte1) & 0x7f) << 1 )  + ((ord(byte2) & 0x80 ) >> 7 )
        fract = ((ord(byte2) & 0x7f) << 16 ) +  (ord(byte3) << 8 ) + ord(byte4)
        if sign == 0: sign_mult = 1.0
        else: sign_mult = -1.0;
        if 0 < expon:
            # note 16777216.0 == 2^24  
            val = sign_mult * (0.5 + (fract/16777216.0)) * pow(2.0, expon - 128.0)   
            f.append(val)
        elif expon == 0 and sign == 0: f.append(0)
        else: f.append(0) # may want to raise an exception here ...
    return f 
def usage():
    lprint ('')
    lprint ('Usage: '+Program_name+' [options] --spec=<spectrofile>')
    lprint ('')
    lprint ('   Available options are:')
    lprint ('       --outdir=<path>   : output directory, if not specified')
    lprint ('                           output goes to current working directory')
    lprint ('       --TE_inc=<number> : TE increment in miliseconds (default '+str(dTE)+'ms)')
    lprint ('                           for SPAR/SDAT input files only')
    lprint ('       --TE1=<number>    : TE1 in miliseconds (default '+str(TE1)+'ms)')
    lprint ('                           twice the intervall between the 90-180 pulses')
    lprint ('                           see --help for details')
    lprint ('       --help (or -h)    : usage and help')
    lprint ('       --version         : version information')
    lprint ('')
def help():
    lprint ('')
    lprint ('the <spectrofile> can be a SPAR/SDAT file') 
    lprint ('   (either one can be specified, but both have to be present)')
    lprint ('or a Philips DICOM spectroscopy file')
    lprint ('   (when exported from the scanner these are called XX*)')
    lprint ('')
    lprint ('TE1 is NOT the first TE of the multi TE acquisition, but')
    lprint ('    the time between the 90\' excitation pulse and the first echo,')
    lprint ('    which is equal to twice the intervall between the 90\' and') 
    lprint ('    the first 180\' pulse of the PRESS pulse sequence.')
    lprint ('    The exact value can only be determined by looking into the timing ')
    lprint ('    details of the pulse sequence')
    lprint ('    As a rule of thumb TE1 = ~0.6*TEmin, where TEmin is the minimum')
    lprint ('    echo time as informed on the info page on the scanner.')
    lprint ('    For TEmin = 35ms --> TE1 = 21ms,') 
    lprint ('    or for protocolls using high gradient settings') 
    lprint ('    TEmin = 30ms --> TE1 = 18ms')
    lprint ('    The TE1 value influences TARQUIN\'s simulation of the') 
    lprint ('    basis functions at different TE\'s, which is important to ')
    lprint ('    correctly model spectral peaks of coupled spins (e.g. Glx) ')
    lprint ('')    
    lprint ('Limitations:')
    lprint (' - SPAR/SDAT files only contain the information of the first echo time')
    lprint ('   per default the software assumes a constant echo time increment of '+str(dTE)+'ms')
    lprint ('   to specify a different value use the --TE_inc option')
    lprint ('   this only applies to SPAR/DAT input')
    lprint ('   for DICOM input TEs are read directly from the DICOM tags')
    lprint (' - to read DICOM format the "pydicom" library is required')
    lprint ('   for installation instructions see http://pydicom.readthedocs.io')     
    lprint (' - to be able to interactively choose the input spectro file')
    lprint ('   the python "tkinter" library is required')
    lprint ('   for installation instructions see http://wiki.python.org/moin/TkInter')    
    lprint ('')

# general initialization stuff   
debug=False; NIFTI_Input=False; SPAR_Input=True
filename=''; workfile=''
slash='/'; 
if sys.platform=="win32": slash='\\' # not really needed, but looks nicer ;)
Program_name = os.path.basename(sys.argv[0]); 
if Program_name.find('.')>0: Program_name = Program_name[:Program_name.find('.')]
basedir = os.getcwd()+slash # current working directory is the default output directory 
for arg in sys.argv[1:]: # look in command line arguments if the output directory specified
    if "--outdir" in arg: basedir = os.path.abspath(arg[arg.find('=')+1:])+slash #
ID = str(random.randrange(1000, 2000));ID=ID[:3] # create 3 digit random ID for logfile 
try: sys.stderr = open(basedir+Program_name+'.log', 'a'); # open logfile to append
except: print('Problem opening logfile: '+basedir+Program_name+'.log'); exit(2)
my_env = os.environ.copy()
FNULL = open(os.devnull, 'w')
# catch signals to be able to cleanup temp files before exit
signal.signal(signal.SIGINT, signal_handler)  # keyboard interrupt
signal.signal(signal.SIGTERM, signal_handler) # kill/shutdown
if  'SIGHUP' in dir(signal): signal.signal(signal.SIGHUP, signal_handler)  # shell exit (linux)
# make tempdir
timestamp=datetime.datetime.now().strftime("%Y%m%d%H%M%S")
tempname='.'+Program_name+'_temp'+timestamp
tempdir=basedir+tempname+slash
if os.path.isdir(tempdir): # this should never happen
    lprint ('ERROR:  Problem creating temp dir (already exists)'); exit(1) 
try: os.mkdir (tempdir)
except: lprint ('ERROR:  Problem creating temp dir: '+tempdir); exit(1) 

# configuration specific initializations
# compare python versions with e.g. if LooseVersion(python_version)>LooseVersion("2.7.6"):
python_version = str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
# sys.platform = [linux2, win32, cygwin, darwin, os2, os2emx, riscos, atheos, freebsd7, freebsd8]
if sys.platform=="win32":
    os.system("title "+Program_name)
    try: resourcedir = sys._MEIPASS+slash # when on PyInstaller 
    except: # in plain python this is where the script was run from
        resourcedir = os.path.abspath(os.path.dirname(sys.argv[0]))+slash; 
    command='attrib'; parameters=' +H "'+tempdir[:len(tempdir)-1]+'"'; 
    run(command, parameters) # hide tempdir
    if pywin32_installed:
        try: # disable console windows close button (substitutes catch shell exit under linux)
            hwnd = win32console.GetConsoleWindow()
            hMenu = win32gui.GetSystemMenu(hwnd, 0)
            win32gui.DeleteMenu(hMenu, win32con.SC_CLOSE, win32con.MF_BYCOMMAND)
        except: pass #silent        
else:
    resourcedir = os.path.abspath(os.path.dirname(sys.argv[0]))+slash;
if TK_installed:        
    TKwindows = tk.Tk(); TKwindows.withdraw() #hiding tkinter window
    TKwindows.update()
    # the following tries to disable showing hidden files/folders under linux
    try: TKwindows.tk.call('tk_getOpenFile', '-foobarz')
    except: pass
    try: TKwindows.tk.call('namespace', 'import', '::tk::dialog::file::')
    except: pass
    try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenBtn', '1')
    except: pass
    try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')
    except: pass
    TKwindows.update()

# parse commandline parameters (if present)
try: opts, args =  getopt( sys.argv[1:],'h',['help','version','spec=','outdir=', 'TE_inc=', 'TE1='])
except:
    error=str(sys.argv[1:]).replace("[","").replace("]","")
    if "-" in str(error) and not "--" in str(error): 
          lprint ('ERROR: Commandline '+str(error)+',   maybe you mean "--"')
    else: lprint ('ERROR: Commandline '+str(error))
    usage(); exit(2)
if len(args)>0: 
    lprint ('ERROR: Commandline option "'+args[0]+'" not recognized')
    lprint ('       (see logfile for details)')
    logwrite ('       Calling parameters: '+str(sys.argv[1:]).replace("[","").replace("]",""))
    usage(); exit(2)  
argDict = dict(opts)
if "--outdir" in argDict and not [True for arg in sys.argv[1:] if "--outdir" in arg]:
    # "--outdir" must be spelled out, getopt also excepts substrings (e.g. "--outd"), but
    # my simple pre-initialization code to get basedir early doesn't
    lprint ('ERROR: Commandline option "--outdir" must be spelled out')
    usage(); exit(2)
if '-h' in argDict: usage(); help(); exit(0)   
if '--help' in argDict: usage(); help(); exit(0)  
if '--version' in argDict: lprint (Program_name+' '+Program_version); exit(0)
if '--spec' in argDict: filename=argDict['--spec']; checkfile(filename)
if '--TE_inc' in argDict: 
    dTE_str=argDict['--TE_inc']
    try: dTE=float(dTE_str)
    except: lprint ('ERROR: problem converting TE_inc to number'); exit(2)
if '--TE1' in argDict: 
    TE1_str=argDict['--TE1']
    try: TE1=float(TE1_str)
    except: lprint ('ERROR: problem converting TE1 to number'); exit(2)
    
#choose file with tkinter
try:
   TKwindows = tk.Tk(); TKwindows.withdraw() #hiding tkinter window
   TKwindows.update()
except: pass
# the following tries to disablle showing hidden files/folders under linux
try: TKwindows.tk.call('tk_getOpenFile', '-foobarz')
except: pass
try: TKwindows.tk.call('namespace', 'import', '::tk::dialog::file::')
except: pass
try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenBtn', '1')
except: pass
try: TKwindows.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')
except: pass
try: TKwindows.update()
except: pass
# this is the actual file open dialogue
Interactive = False
if TK_installed:
   # Choose Spectro file
    if filename == "": # use interactive input if not specified in commandline
        filename = askopenfilename(title="Choose Spectro file")
        if filename == "": lprint ('ERROR:  No Spectro input file specified'); exit(2)
        Interactive = True
    TKwindows.update()
else:
    if filename == "": 
        lprint ('ERROR:  No Spectro input file specified')
        lprint ('        to interactively choose input files you need tkinter')
        lprint ('        on Linux try "yum install tkinter"')
        lprint ('        on MacOS install ActiveTcl from:')
        lprint ('        http://www.activestate.com/activetcl/downloads')  
        usage()
        exit(2)
filename = os.path.abspath(filename)


# ----- start to really do something -----
lprint ('Starting '+Program_name+' '+Program_version)
logwrite ('Calling sequence    '+' '.join(sys.argv))
logwrite ('OS & Python version '+sys.platform+' '+python_version)
logwrite ('tkinter & pydicom   '+str(TK_installed)+' '+str(pydicom_installed))

# read data
if isDICOM (filename):  # read DICOM  
    try: Dset = dicom.read_file(filename)
    except: lprint ('ERROR: reading DICOM file'); exit(1)
    # do some checks
    try: Modality=str(Dset.Modality) # must be MR           
    except: lprint ('ERROR: Unable to determine DICOM Modality'); exit(1)
    if Modality!='MR': lprint ('DICOM Modality not MR'); exit(1)
    try: Manufacturer=str(Dset.Manufacturer) # currently Philips only
    except: lprint ('ERROR: Unable to determine Manufacturer'); exit(1)
    if not Manufacturer.find("Philips")>=0: 
        lprint ('ERROR: Currently only Philips DICOM implemented'); exit(1)
    try: ImageType=str(Dset.ImageType) # sanity check: spectroscopy   
    except: 
        lprint ('ERROR: Unable to determine if DICOM containes spectroscopy data'); exit(1)
    if not ImageType.find("SPECTROSCOPY")>=0: 
        lprint ('ERROR: DICOM file does not contain spectroscopy data'); exit(1)
    # start reading data from DICOM
    ReIm=2 # real and imagininary parts
    try: samples=Dset.SpectroscopyAcquisitionDataColumns # n_points in time/frequency domain
    except: lprint ('ERROR: reading number of samples from DICOM file'); exit(1)
    try: rows=Dset[0x2001,0x1081].value #NumberOfDynamicScans
    except: lprint ('ERROR: reading number of dynamics from DICOM file'); exit(1)
    if rows<=1: lprint ('ERROR: not a multi TE aquisition'); exit(1)
    if rows<=min_nTE: lprint ('ERROR: minimum number of'+str(min_nTE)+'TEs required'); exit(1)
    # check if number of points is correct
    ndata_points=len (numpy.asarray(Dset[0x5600,0x0020].value))
    if ndata_points==2*rows*samples*ReIm: 
        ActRef=2 # two spectra, actual and water, this is the normal case
    elif ndata_points==rows*samples*ReIm: 
        ActRef=1 # only one spectrum
    else: 
        lprint ('ERROR: Unexpected number of total datapoints'); exit(1)  
    # read TE values
    TE = numpy.zeros(shape=(rows))
    for i in range (0, rows): 
        TE[i]=Dset.PerFrameFunctionalGroupsSequence[i].MREchoSequence[0].EffectiveEchoTime
    if TE[1:].tolist()==TE[:-1].tolist(): # all elements equal (to check TE=numpy.array([1,1,1]))
        lprint ('ERROR:: all TEs are equal'); exit(1)
    #read data
    spectro_rawdata = numpy.asarray(Dset[0x5600,0x0020].value)
else: # read SDAT
    # find SPAR/SDAT pair
    path=os.path.dirname(filename)
    name=os.path.splitext(os.path.basename(filename))[0]
    ext=os.path.splitext(os.path.basename(filename))[1]
    if ext.lower() == ".spar": 
      SPARfile=filename 
      SDATfile=[f for f in os.listdir(path) if f.lower().endswith('.sdat') and f.startswith(name)]
      SDATfile=SDATfile[0] # may want to raise an exception here ...
    elif ext.lower() == ".sdat":
      SDATfile=filename 
      SPARfile=[f for f in os.listdir(path) if f.lower().endswith('.spar') and f.startswith(name)]
      SPARfile=SPARfile[0] # may want to raise an exception here ...
    else: lprint ('ERROR: file extension should be SDAT/SPAR'); exit(1)
    # open SPAR
    try: input = open(SPARfile, "r").readlines()
    except: lprint ('ERROR: reading SPAR file'); exit(1)
    TE_first = float(_get_from_SPAR (input, 'echo_time'))
    samples = int(_get_from_SPAR (input, 'samples'))
    rows = int(_get_from_SPAR (input, 'rows'))
    ActRef = int(_get_from_SPAR (input, 'mix_number'))
    if ActRef != 1: 
       lprint ('ERROR: SPAR/SDAT file seems to be a reference spectrum, choose an actual spectrum'); 
       exit(1)
    # guess missing data not contained in the SPAR file
    TE=numpy.linspace(TE_first,TE_first+rows*dTE,rows,endpoint=False) # SPAR contains only first TE
    ReIm=2   # real and imagininary parts
logwrite ('TE1 = '+str(TE1))
logwrite ('Reading File '+filename)   
lprint ('') # spacer
    
# DICOM bruteforce remove REF spectrum
if isDICOM(filename):
    if ActRef==2:
        spectro_rawdata = numpy.asarray(Dset[0x5600,0x0020].value)
        spectro_rawdata = spectro_rawdata [0:spectro_rawdata.size/2] # remove REF
        Dset[0x5600,0x0020].value = spectro_rawdata.tolist()
        
        
# start processing with TARQUIN
##delete (tempdir+'tarquin_T2fit.csv') # cleanup
##delete (tempdir+'tarquin_T2fit.txt') # cleanup
##delete (tempdir+'avlist.csv')        # cleanup
##delete (tempdir+workfile)            # cleanup
signal=numpy.zeros((len(metabolites),rows))
std_dev=numpy.zeros((len(metabolites),rows))
for n_TE in range(rows):
   space=''
   if n_TE<9: space=' '
   lprint ('Processing spectrum '+space+str(n_TE+1)+' of '+str(rows))
   avfile = open(tempdir+'avlist.csv', 'w')
   avfile.write(str(n_TE+1))  
   avfile.close()
   if isDICOM(filename): 
      workfile='Spectro_noREF.DCM'
      try: Dset.PerFrameFunctionalGroupsSequence[0].MREchoSequence[0].EffectiveEchoTime = TE[n_TE]
      except: lprint ("ERROR: problem modifying TE"); exit(1)
      Dset.save_as(tempdir+workfile) # write dicom data back to file
      arguments=' --input "'+tempdir+workfile+'"'
   else:
      arguments=' --input "'+filename+'"'
      arguments+=' --format philips'
      arguments+=' --echo '+str(TE[n_TE]/1000.) # seems to have no effect on DICOM files
   arguments+=' --av_list "'+tempdir+'avlist.csv"' 
   arguments+=' --output_csv  "'+tempdir+'tarquin_T2fit.csv"'
   arguments+=' --output_txt  "'+tempdir+'tarquin_T2fit.txt"'
   arguments+=' --ref 4.66 --max_metab_shift 0.015'
   arguments+=' --auto_phase true --dyn_freq_corr true'
   arguments+=' --te1 '+str(TE1/1000.)
   arguments+=' --start_pnt 20 --ref_signals 1h_naa --dref_signals 1h_naa --pul_seq press --int_basis 1h_brain_full'
   run (resourcedir+'tarquin', arguments)
   with open(tempdir+'tarquin_T2fit.csv', 'r') as csvfile:
      data_iter = csv.reader(csvfile, delimiter=',')
      data = [data for data in data_iter]
      # read header
      header = data[1]
      metabolite_idx = []
      for i in range(len(metabolites)):
         try: metabolite_idx.append (header.index(metabolites[i]))
         except:
            lprint ('ERROR: metabolite '+metabolites[i]+' not found in tarquin output')
            exit(1)
      # read signal & stddev
      for i in range(len(metabolites)):
         data_row = data[2]
         signal[i,n_TE] = data_row[metabolite_idx[i]]
         data_row = data[5]
         std_dev[i,n_TE] = data_row[metabolite_idx[i]]
         # use relative error, otherwise use curve_fit (...., absolute_sigma=True)
         std_dev[i,n_TE] /= signal[i,n_TE] 
   
#fit & write results
lprint ('') # spacer
stp=''; space = ' ' # for name collision detection
if os.path.isfile(basedir+Program_name+'_Results.txt'): stp=timestamp+ID+'_'
f = open(basedir+stp+Program_name+'_Results.txt', 'w')
f.write(Program_name+space+Program_version+' Results:\n')
for i in range(len(metabolites)): f.write(metabolites[i]+' \t')
f.write('Spectro_File \tbasedir\n')
# start fitting
for i in range(len(metabolites)):
   initial_conditions=[max(signal[i,:])/2, 100]
   pars,covar = curve_fit(expdecay,TE,signal[i,:], p0=initial_conditions,  
                maxfev=100*(len(TE)+1),sigma=std_dev[i,:]) 
   T2=pars[1] 
   dT2=covar[1,1]
   chi2 = sum((signal[i,:]-expdecay(TE,*pars))**2/expdecay(TE,*pars))/(len(TE) - 2)
   T2_str  = "%.2f" % T2
   dT2_str = "%.2f" % dT2
   metabolite_name=metabolites[i]
   if metabolite_name[0]=='T': metabolite_name=metabolite_name[1:]
   while len(metabolite_name)<3: metabolite_name += ' '
   lprint (metabolite_name+':  T2='+T2_str)
   logwrite ('StdDevT2='+dT2_str)
   logwrite ('chi2    ='+str(chi2))
   logwrite (' \t TE \t'+'Signal('+metabolite_name+') \t'+'StdDev('+metabolite_name+')')
   f.write(T2_str+' \t')
   for j in range(rows):
      logwrite (' \t'+str(TE[j])+' \t'+str(signal[i,j])+' \t'+str(std_dev[i,j]))
   logwrite ('') # spacer
f.write(filename+' \t')
f.write(basedir+' \n') 
f.close()   
lprint ('') # spacer

#delete tempdir
try: shutil.rmtree(tempdir)
except: pass # silent
lprint ('done\n')
sys.stderr.close() # close logfile

#reenable console windows close button
if pywin32_installed:
    try:
        hwnd = win32console.GetConsoleWindow()
        hMenu = win32gui.GetSystemMenu(hwnd, 1)
        win32gui.DeleteMenu(hMenu, win32con.SC_CLOSE, win32con.MF_BYCOMMAND)
    except: pass #silent

#pause
if Interactive:
    if sys.platform=="win32": os.system("pause") # windows
    else: 
        #os.system('read -s -n 1 -p "Press any key to continue...\n"')
        import termios
        print("Press any key to continue...")
        fd = sys.stdin.fileno()
        oldterm = termios.tcgetattr(fd)
        newattr = termios.tcgetattr(fd)
        newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
        termios.tcsetattr(fd, termios.TCSANOW, newattr)
        try: result = sys.stdin.read(1)
        except IOError: pass
        finally: termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)   


