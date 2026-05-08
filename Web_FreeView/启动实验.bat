@echo off
setlocal
cd /d %~dp0

set PORT=8000
set URL=http://127.0.0.1:%PORT%/
set PYTHON_EXE=

start "" "%URL%"

if exist ".\bin\python\python.exe" (
	echo [INFO] using bundled runtime .\bin\python\python.exe
	set PYTHON_EXE=.\bin\python\python.exe
)

if not defined PYTHON_EXE (
	where python >nul 2>nul
	if %errorlevel%==0 (
		echo [INFO] using system python
		set PYTHON_EXE=python
	)
)

if not defined PYTHON_EXE (
	echo [ERROR] python not found. Please place a Python runtime in .\bin\python or install Python.
	pause
	goto :end
)

echo [INFO] starting local server on port %PORT%
%PYTHON_EXE% -u ".\bin\start_server.py" --port %PORT%

:end
endlocal