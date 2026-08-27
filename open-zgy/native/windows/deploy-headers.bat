if [%~1]==[] goto usage

set OUTDIR=%1
xcopy /i /y ..\src\*.h %OUTDIR%include\openzgy
exit /b 0

:usage
@echo Usage: %0 ^<OUTDIR^>
@exit /b 1
