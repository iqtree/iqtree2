# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# For each target create a dummy rule so the target does not have to exist
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/libz.dylib:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/libz.dylib:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/libz.dylib:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/libz.dylib:


# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.example.Debug:
PostBuild.zlib.Debug: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/example
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/example:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/example


PostBuild.minigzip.Debug:
PostBuild.zlib.Debug: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/minigzip
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/minigzip:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/minigzip


PostBuild.zlib.Debug:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/libz.dylib:
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Debug/libz.dylib


PostBuild.zlibstatic.Debug:
PostBuild.example.Release:
PostBuild.zlib.Release: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/example
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/example:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/example


PostBuild.minigzip.Release:
PostBuild.zlib.Release: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/minigzip
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/minigzip:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/minigzip


PostBuild.zlib.Release:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/libz.dylib:
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/Release/libz.dylib


PostBuild.zlibstatic.Release:
PostBuild.example.MinSizeRel:
PostBuild.zlib.MinSizeRel: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/example
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/example:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/example


PostBuild.minigzip.MinSizeRel:
PostBuild.zlib.MinSizeRel: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/minigzip
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/minigzip:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/minigzip


PostBuild.zlib.MinSizeRel:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/libz.dylib:
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/MinSizeRel/libz.dylib


PostBuild.zlibstatic.MinSizeRel:
PostBuild.example.RelWithDebInfo:
PostBuild.zlib.RelWithDebInfo: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/example
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/example:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/example


PostBuild.minigzip.RelWithDebInfo:
PostBuild.zlib.RelWithDebInfo: /Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/minigzip
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/minigzip:\
	/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/libz.dylib
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/minigzip


PostBuild.zlib.RelWithDebInfo:
/Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/libz.dylib:
	/bin/rm -f /Users/minh/projects/iqtree-cmake/zlib-1.2.7/RelWithDebInfo/libz.dylib


PostBuild.zlibstatic.RelWithDebInfo:
