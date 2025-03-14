:: @echo off
setlocal enabledelayedexpansion

:: Set the output library name
set OUTPUT_LIB=libiqtree2full.a
set OUTPUT_LIB2=libiqtree2full.lib

:: Define the temporary folder
set TEMP_DIR=temp_objs

:: Cleanup any existing temporary directory
if exist %TEMP_DIR% rmdir /s /q %TEMP_DIR%
mkdir %TEMP_DIR%

:: Iterate through all .a files passed as arguments
for %%L in (%*) do (
    set LIB_NAME=%%~nL
    set EXT=%%~xL
    mkdir %TEMP_DIR%\!LIB_NAME!

    if "!EXT!"==".dll" (
        echo Creating import library from !LIB_NAME!.dll...

        :: Step 1: Generate a .def file
        gendef "%%L"

        :: Step 2: Use dlltool to create an import library
        dlltool -d "!LIB_NAME!.def" -l "!LIB_NAME!.a" -k

        :: Step 3: Extract object files from the new .a file
        ar x "!LIB_NAME!.a"
        move *.o %TEMP_DIR%\!LIB_NAME!\ >nul 2>&1
        move *.obj %TEMP_DIR%\!LIB_NAME!\ >nul 2>&1

    ) else if "!EXT!"==".a" (
        echo Extracting object files from !LIB_NAME!.a...
        ar x "%%L" >nul 2>&1
        move *.o %TEMP_DIR%\!LIB_NAME!\ >nul 2>&1
        move *.obj %TEMP_DIR%\!LIB_NAME!\ >nul 2>&1

    ) else (
        echo Skipping unsupported file type: %%L
    )
)

:: Merge all extracted .obj files into the final static library
cd %TEMP_DIR%

(for /r %%F in (*.o *.obj) do (
    set "p=%%F"
    set "p=!p:\=/!"
    echo !p!
)) > filelist.txt

ar rcs ../%OUTPUT_LIB% @filelist.txt
lib /out:../%OUTPUT_LIB2% @filelist.txt

cd ..

:: Cleanup extracted files
rmdir /s /q %TEMP_DIR%