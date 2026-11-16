@ECHO OFF
setlocal

pushd %~dp0

REM Command file for Sphinx documentation

if not defined SPHINXBUILD (
	set SPHINXBUILD=python -msphinx
)
set SOURCEDIR=.
set BUILDDIR=../docs
set SPHINXPROJ=schimpy

if "%1" == "" goto help

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The Sphinx module was not found. Make sure you have Sphinx installed,
	echo.then set the SPHINXBUILD environment variable to point to the full
	echo.path of the 'sphinx-build' executable. Alternatively you may add the
	echo.Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	set BUILD_STATUS=1
	goto end
)

%SPHINXBUILD% -M %1 "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS%
set BUILD_STATUS=%errorlevel%
goto end

:help
%SPHINXBUILD% -M help "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS%
set BUILD_STATUS=%errorlevel%

:end
popd
exit /b %BUILD_STATUS%
