@echo off
REM Build lock — only one build at a time across all scripts
set LOCKFILE=F:\Claude\.build_lock
if exist "%LOCKFILE%" (
    echo ERROR: Another build is already running. Lockfile: %LOCKFILE%
    exit /b 1
)
echo %DATE% %TIME% > "%LOCKFILE%"

REM Put Strawberry Perl first on PATH (OpenSSL needs Windows-native perl)
set PATH=F:\Claude\strawberry-perl\perl\bin;F:\Claude\strawberry-perl\c\bin;%PATH%
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
set CMAKE_POLICY_VERSION_MINIMUM=3.5
cd /d F:\Claude\OrcaSlicer\deps
if not exist build mkdir build
cd build
cmake ../ -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target deps -j 12

del "%LOCKFILE%" 2>NUL
