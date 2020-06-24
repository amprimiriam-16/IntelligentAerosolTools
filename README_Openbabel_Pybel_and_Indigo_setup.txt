0.) Note that the following steps are outlined for a Windows 64bit installation; similar steps may need to be taken on a Linux machine.

1.) Openbabel and pybel installation:
------------------------------------
- first make sure that the python 'pip' is installed and up to date:
	in a command prompt with administrator rights (click Run as Administrator) type:  python -m pip install --upgrade pip

- install openbabel 3.0.0 (GUI) 64bit version for Windows (download the installer file)

- install python 3.7.7 64bit version (somehow the python bindings of openbabel don't work yet with python 3.8+)

- install python bindings: run in a command prompt with administrator rights:  pip install -U openbabel
	
	and follow other setup/testing advice here:  https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html
	
- commands like obabel and 'python form openbabel import pybel' should work without error now.


2.) Indigo installation
-----------------------
For use of the Indigo toolkit, the indigo package and python bindings need to be installed:
- run in a command prompt with administrator rights:  
	pip install epam.Indigo
	
	


 