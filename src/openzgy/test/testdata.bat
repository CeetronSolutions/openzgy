if [%~1]==[] goto usage

set OUTDIR=%1

7z x -y -o%OUTDIR% ..\..\testdata\Empty-v1.zgy.bz2
7z x -y -o%OUTDIR% ..\..\testdata\Empty-v3.zgy.bz2
7z x -y -o%OUTDIR% ..\..\testdata\Fancy-int8.zgy.bz2
7z x -y -o%OUTDIR% ..\..\testdata\EmptyOldFile.zgy.bz2

copy /B %OUTDIR%Empty-v1.zgy+,, %OUTDIR%Empty-v1.zgy
copy /B %OUTDIR%Empty-v3.zgy+,, %OUTDIR%Empty-v3.zgy
copy /B %OUTDIR%Fancy-int8.zgy+,, %OUTDIR%Fancy-int8.zgy
copy /B %OUTDIR%EmptyOldFile.zgy+,, %OUTDIR%EmptyOldFile.zgy

exit /b 0

:usage
@echo Usage: %0 ^<OUTDIR^>
@exit /b 1
