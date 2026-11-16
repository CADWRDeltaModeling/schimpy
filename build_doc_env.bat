@echo off
setlocal
if "%~1"=="--help" goto help
if "%~1"=="-h" goto help
if not "%~1"=="" (
    echo Usage: build_doc_env.bat
    exit /b 2
)

pushd "%~dp0"
call conda env create -n schimpy_doc -f environment.yml
set RESULT=%errorlevel%
if not "%RESULT%"=="0" goto end
call conda env update -n schimpy_doc -f environment-doc.yml
set RESULT=%errorlevel%
if not "%RESULT%"=="0" goto end
call conda activate schimpy_doc
set RESULT=%errorlevel%
if not "%RESULT%"=="0" goto end
python -m pip install --no-deps -e .
set RESULT=%errorlevel%
if not "%RESULT%"=="0" goto end
echo Documentation environment ready.
echo Next: conda activate schimpy_doc
echo Then: cd docsrc
echo       make html

:end
popd
exit /b %RESULT%

:help
echo Usage: build_doc_env.bat
echo Create schimpy_doc from environment.yml plus environment-doc.yml,
echo then install this checkout editable. Run in a Conda-initialized terminal.
echo An existing environment is not deleted or replaced.
exit /b 0
