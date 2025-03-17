:: @echo off
setlocal enabledelayedexpansion

:: Build the DLL file and the DEF file
g++ -shared -o libiqtree2.dll %1 -Wl,--output-def,libiqtree2.def -liqtree2 -L. -lws2_32 -lz -lpthread

:: Create the LIB file
lib /def:libiqtree2.def
