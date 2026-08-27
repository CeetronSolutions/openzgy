if [%~1]==[] goto usage
if [%~2]==[] goto usage
if [%~3]==[] goto usage
if [%~4]==[] goto usage

set OUTDIR=%1
set PLATFORM=%2
set CONFIG=%3
set MYTEMP=%4%zfpbuildnative
@rem set CMAKE="c:\Program Files\CMake\bin\cmake.exe"
set CMAKE=cmake

@rem ... Choosing to store the ZFP source code as a single file in git to make
@rem ... sure we use the exact version from the repos.
rmdir /s /q %MYTEMP%
7z x -o%MYTEMP% ..\..\external\zfp-0.5.5.tar.gz
7z x -o%MYTEMP% %MYTEMP%\zfp-0.5.5.tar
md %MYTEMP%\build
cd %MYTEMP%\build

@rem ... I am unsure whether cmake picks up the right compiler based on the
@rem ... environment. I might need to be explicit. Also note that -A is Note
@rem ... valid for VS2015 Win64, as architecture is already provided.
@rem ... So this won't work for win32.
%CMAKE% -G "Visual Studio 17 2022" %MYTEMP%\zfp-0.5.5
@rem %CMAKE% -A %PLATFORM% %MYTEMP%\zfp-0.5.5
%CMAKE% --build . --config %CONFIG%
copy bin\%CONFIG%\zfp.dll %OUTDIR%
copy bin\%CONFIG%\zfp.pdb %OUTDIR%
copy lib\%CONFIG%\zfp.lib %OUTDIR%

@rem ... Note that CMake has been better integrated in newer versions of
@rem ... visual studio. Consider leveraging that instead of this script.
@rem ... Also, Azure DevOps also has some CMake integration but if I use
@rem ... that I might make it harder to build locally.
@rem ... https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2015
@rem ... https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2019

@rem ... How it works today:
@rem ... I have added  a separate "ZFP" project type "Utility" to the solution.
@rem ... Solution->Project Dependencies... (make OpenZGY depend on ZFP)
@rem ... build-zfp.bat->Properties "Custom Build Tool" run this script
@rem ... OpenZGY->Properties->Linker->Input->AdditionalDependencies (zfp.lib)
@rem ... OpenZGY->src\impl\compress.zfp->Properties->
@rem ...     C/C++->Additional Include Directories->$(???)zfp-0.5.5\include
 
 exit /b 0

:usage
@echo Usage: %0 ^<OUTDIR^> ^<PLATFORM^> ^<CONFIG^> ^<TEMPDIR^>
@exit /b 1
