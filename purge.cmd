@echo off
for /d /r %%G in ("slprj") do rd /s /q "%%~G"
del /s /q *.slxc
del /s /q *.slx.autosave