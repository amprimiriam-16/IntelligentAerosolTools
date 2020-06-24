@ECHO OFF
SETLOCAL

ECHO ----

SET _inputPath=./InputFiles/

:: set variables for this calculation case (input file name):
SET _fileN=smiles_1999.txt

:: run the python code with the chosen input properties:
python .\SMILES_to_AIOMFAC_input.py  %_inputPath%%_fileN%

ECHO -------------------------------------
ECHO The python batch job is finished. 
ECHO -------------------------------------

pause
