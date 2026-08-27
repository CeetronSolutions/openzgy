
@rem ... Usage: test-wrapper.bat WRAPPER_DIR WHEEL HAVE_SD

@echo RUN %0 %1 %2 %3
@echo AZURE_BUILDID=%AZURE_BUILDID%
@echo BUILD_BUILDID=%BUILD_BUILDID%

@echo SD_TARGET_AUDIENCE: '%SD_TARGET_AUDIENCE%'
@echo SDAPI_SAUTH_TOKEN_SERVICE_URL: '%SDAPI_SAUTH_TOKEN_SERVICE_URL%'

echo off
set SD_OK=maybe
if [%~3] == [yes] (
  echo SD tests were requested.
  if not [%OPENZGY_SDURL%] == [] (
    echo OPENZGY_SDURL is set
    if not [%OPENZGY_SDAPIKEY%] == [] (
      echo OPENZGY_SDAPIKEY is set
      if not [%OPENZGY_TOKEN%] == [] (
        echo OPENZGY_TOKEN is set
        set SD_OK=yes
      )
    )
  )
)
echo on

if [%~3] == [yes] (
  if [%SD_OK%] == [yes] (
    echo SD tests will be run.
  ) else (
    echo WARNING: SD tests were requested but secrets are missing.
    echo SD_OK=%SD_OK%.
    set SD_OK=no
  )
) else (
  echo WARNING: SD tests are disabled.
  echo SD_OK=%SD_OK%.
)

if [%AZURE_BUILDID%] == [] set AZURE_BUILDID=%BUILD_BUILDID%

set WRAPPER_SRC=%~f1
set WRAPPER_WHEEL=%~f2

@echo Wrapper: setting up test environment
cd %WRAPPER_SRC%\test
rd /s /q venv
python -m venv venv
call venv\Scripts\activate.bat
python -m pip install  --upgrade pip
python -m pip install  --upgrade numpy
@echo Wrapper: installing %WRAPPER_WHEEL%
python -m pip install  --upgrade %WRAPPER_WHEEL%
python -m pip list

@echo OPENZGY_SDK=%OPENZGY_SDK%
@rem ... need backslashes ... dir %OPENZGY_SDK%
@echo OpenZGY site package
dir venv\lib\site-packages\openzgycpp
@echo Installed wheel
dir %WRAPPER_WHEEL%

@echo Wrapper: running test.py
python test.py
if ERRORLEVEL 1 echo Wrapper: test.py completed with errorlevel %ERRORLEVEL%.
if ERRORLEVEL 1 exit /b %ERRORLEVEL%

@echo Wrapper: running test_pickle.py
python test_pickle.py
if ERRORLEVEL 1 echo Wrapper: test_pickle.py completed with errorlevel %ERRORLEVEL%.
if ERRORLEVEL 1 exit /b %ERRORLEVEL%

@echo Wrapper: running test_black.py local
python test_black.py local
if ERRORLEVEL 1 echo Wrapper: test_black.py local completed with errorlevel %ERRORLEVEL%.
if ERRORLEVEL 1 exit /b %ERRORLEVEL%

@echo Wrapper: maybe do the SD tests also (%SD_OK%) ?
if [%SD_OK%] == [yes] (
  @echo Wrapper: running test_black.py sd
  python test_black.py sd
  if ERRORLEVEL 1 echo Wrapper: test_black.py sd completed with errorlevel %ERRORLEVEL%.
  if ERRORLEVEL 1 exit /b %ERRORLEVEL%
) else (
  @echo Wrapper: skipping test_black.py sd
)
@echo Wrapper: DONE TESTING - ALL OK.

exit /b 0
