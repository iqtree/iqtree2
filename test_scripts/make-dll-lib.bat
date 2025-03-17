:: @echo off
setlocal enabledelayedexpansion

:: Build the DLL file and the DEF file
g++ -shared -o iqtree2.dll %1 -Wl,--output-def,iqtree2.def -liqtree2 -L. -lws2_32 -lz -lpthread

:: Create the LIB file
lib /def:iqtree2.def
