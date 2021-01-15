::to use, run "" .\upload.bat d### "" in the terminal, where d### is the distribution
::automatically creates wheel and tar, uploads to pypi and anaconda. 
@ECHO OFF



:: set argument equal to variable, create the path to distribution, create the wheel and tar
set dist_folder=%1
set "distpath=dist/%dist_folder%"
py setup.py sdist bdist_wheel --dist-dir %distpath%


:: UPLOAD TO PYPI
set "distpath=%distpath%/*"
py -m twine upload --repository pypi %distpath%


:: STRING HANDLING AND FILENAME CREATION FOR ANACONDA UPLOAD
set dist=%dist_folder:~1%
::this makes d### -> ###
set first=%dist:~0,1%
set second=%dist:~1,1%
set third=%dist:~2,2%
set letter=%dist:~-1%
IF %letter% == a set "third=%third%0"
IF %letter% == b set "third=%third%0"
IF %letter% == c set "third=%third%0"
set "conda_dist=gas_dynamics-%first%.%second%.%third%.tar.gz"


:: ACTIVATE CONDA AND UPLOAD TO CONDA
chdir dist
call C:\ProgramData\Anaconda3\Scripts\activate.bat
call anaconda upload %conda_dist%
call conda deactivate