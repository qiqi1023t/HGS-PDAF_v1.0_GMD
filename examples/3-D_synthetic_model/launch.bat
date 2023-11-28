@echo off

python initialize.py

for /l %%i in (1,1,3) do (
cd .\worker_%%i/
start launch.bat
cd ..
)
