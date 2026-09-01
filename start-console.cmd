@echo off
rem Double-click this file to start the Fermentation Console.
rem
rem It does three things: check that uv is installed, let uv install what the console needs,
rem and hand over to app/start_console.py, which picks a free port and opens the browser.
rem Nothing here installs anything system-wide, and nothing is uploaded anywhere.

cd /d "%~dp0"
title Fermentation Console

where uv >nul 2>nul
if errorlevel 1 (
    echo.
    echo   The console is started by a tool called uv, which is not installed on this machine.
    echo.
    echo   Install it once from:  https://docs.astral.sh/uv/getting-started/installation/
    echo   Then double-click this file again.
    echo.
    pause
    exit /b 1
)

echo.
echo   Getting the console ready. The first time takes a minute or two while the
echo   pieces are downloaded; after that it is a few seconds.
echo.

uv run --group ui python app/start_console.py
set EXITCODE=%ERRORLEVEL%

if not "%EXITCODE%"=="0" (
    echo.
    echo   The console stopped with a problem. The lines above say what it was.
    echo.
)
pause
exit /b %EXITCODE%
