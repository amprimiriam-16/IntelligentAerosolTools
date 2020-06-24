@ECHO OFF
SETLOCAL

ECHO ----

SET _inputPath=./InputFiles/

:: set variables for this calculation case (input file name and temperature):
SET _file=smiles_1999.txt
SET _TempK=298.15

:: run the python code with the chosen input properties:
python .\PureCompProperties_for_LumpingSchemes.py  %_inputPath%%_file%  %_TempK%

ECHO 
ECHO -------------------------------------
ECHO The python batch job is finished. 
ECHO -------------------------------------
ECHO.

pause