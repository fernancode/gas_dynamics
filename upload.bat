@ECHO OFF
::to use, run "" .\upload.bat d### "" in the terminal

set dist_folder=%1
:: set dist = first input from command line


:: CREATE THE WHEEL AND TAR.GZ
set "distpath=dist/%dist_folder%"
py setup.py sdist bdist_wheel --dist-dir %distpath%


:: UPLOAD TO PYPI
set "distpath=%distpath%/*"
py -m twine upload --repository pypi %distpath%


:: UPLOAD TO PYPI
:: STRING HANDLING AND FILENAME CREATION
set dist=%dist_folder:~1%
::this makes d### -> ###
:: CREATE THE FILENAME FOR THE CONDA UPLOAD SOME STRING HANDLING
set first=%dist:~0,1%
set second=%dist:~1,1%
set third=%dist:~2,2%
set letter=%dist:~-1%
IF %letter% == a set "third=%third%0"
IF %letter% == b set "third=%third%0"
IF %letter% == c set "third=%third%0"
set "conda_dist=gas_dynamics-%first%.%second%.%third%.tar.gz"
C:
:: ACTIVATE CONDA AND UPLOAD TO CONDA
chdir C:\ProgramData\Microsoft\Windows\"Start Menu"\Programs\"Anaconda3 (64-bit)"\
::call "Anaconda Prompt (Anaconda3).lnk"
::call D:
::call chdir D:\Documents\Python\gas_dynamics\
::call anaconda upload dist/conda_dist
::call deactivate