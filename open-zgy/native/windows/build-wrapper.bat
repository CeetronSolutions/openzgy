
@rem ... Usage: build-wrapper.bat WRAPPER_DIR $(OutputPath) $(Platform) $(Config)

@echo RUN %0 %1 %2 %3 %4
@echo AZURE_BUILDID=%AZURE_BUILDID%
@echo BUILD_BUILDID=%BUILD_BUILDID%
if [%AZURE_BUILDID%] == [] set AZURE_BUILDID=%BUILD_BUILDID%

python --version
python -c "import sys;print(sys.executable, *sys.path, sep='\n')"

@rem ... Trying to work around spurious failures (cannot find python3.dll)
@rem ... when building on some random developer's machine. On Azure DevOps
@rem ... UsePythonVersion sets this variable and does other stuff to make
@rem ... Python work. Not sure whether it helps to set it here. And
@rem ... definitely don't check in this code. See also the kludge in
@rem ... wrapper/setup.py
@rem set USEPYTHONVERSION_PYTHONLOCATION=C:\hostedtoolcache\windows\Python\3.6.8\x64
@rem echo ENV
@rem set
@rem echo @END

cd %1

if [%~3] == [x64] (
    goto platform_ok
)
(
    @echo Platform %3 is not supported.
    exit /b 1
)
:platform_ok

if [%~4] == [Release] (
    echo [build_ext] > setup.cfg
    goto config_ok
)
if [%~4] == [Debug] (
    echo [build_ext] > setup.cfg
    echo debug = 1  >> setup.cfg
    goto config_ok
)
(
    @echo Config %3 is not supported.
    exit /b 1
)
:config_ok
type setup.cfg

rd /s /q test-venv
rd /s /q build
rd /s /q dist
rd /s /q OpenZgyBindings.egg-info

@rem setup links to make the source tree look ok to setuptools.
del openzgy-sdk\__init__.py
rd  openzgy-sdk\include
rd  openzgy-sdk
mklink /j openzgy-sdk ..\build\deploy\native\%3\%4
mklink /j openzgy-sdk\include ..\native\src
copy openzgycpp\__init__.py openzgy-sdk\__init__.py

python -m venv test-venv
call test-venv\Scripts\activate.bat
python -m pip install  --upgrade pip
python -m pip install  --upgrade setuptools wheel

@echo EXEC: python setup.py --no-user-cfg sdist bdist_wheel
python setup.py --no-user-cfg sdist bdist_wheel
dir dist
dir build
@echo EXEC: copy /b /y dist\* %2
md %2
copy /b /y dist\* %2

@rem force rebuild when doing the generic package.
rd /s /q build
rd /s /q dist
rd /s /q OpenZgyBindings.egg-info

@echo EXEC: (vanilla) python setup.py --no-user-cfg sdist bdist_wheel
set OPENZGY_VANILLA=1
python setup.py --no-user-cfg sdist bdist_wheel
set OPENZGY_VANILLA=
dir dist
dir build
@echo EXEC: copy /b /y dist\* %2\..\..\generic\%4
md %2\..\..\generic\%4
copy /b /y dist\* %2\..\..\generic\%4

@rem clean up the links
del openzgy-sdk\__init__.py
rd  openzgy-sdk\include
rd  openzgy-sdk

echo > %2\.\OpenZgyBindings.timestamp
@echo DONE WITH PYTHON

rd /s /q test-venv
rd /s /q build
rd /s /q dist
rd /s /q OpenZgyBindings.egg-info
exit /b 0
