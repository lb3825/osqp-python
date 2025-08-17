@echo off
setlocal enabledelayedexpansion

set PROBNAME=%1
set STARTDIR=%CD%

cd examples
@REM py maros_meszaros_benchmarks.py !PROBNAME!
cd ..\..

cd osqp
py .\plot\plot.py !PROBNAME! 3

cd /d !STARTDIR!