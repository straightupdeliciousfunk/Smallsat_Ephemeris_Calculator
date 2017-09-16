@echo off
:loop
python C:\Users\Tim\workspace\ephemeris_calculator\ephemeris.py
timeout /t 60 >nul
goto loop
