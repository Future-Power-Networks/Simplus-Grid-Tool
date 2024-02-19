@echo off
cd ..
for /d /r %%G in ("slprj","*_grt_rtw") do rd /s /q "%%~G"
del /s /q *.slxc
del /s /q *.slx.autosave