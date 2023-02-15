@echo off

:: Run this script to compile!

:: But first, make sure to substitute this path for wherever your vcvarsall.bat is, which will depend on the version of your MSVC compiler.
call "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" x64

IF NOT EXIST .\build mkdir .\build
pushd .\build

cl ..\main.cpp -Femain.exe -Z7 -nologo -link User32.lib Opengl32.lib Gdi32.lib

popd

pause
