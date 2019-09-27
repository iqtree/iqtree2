#! /bin/bash

# This script is for when you intend to actually work on the code.
# Do not use it, if all you want to do is build the library for
# actual use. If that is what you want, just go through the usual
# way steps for building a cmake-project.

ln -s "../../scripts/pre-commit" ".git/hooks/pre-commit"

if command -v "ninja" >/dev/null; then
	echo "found ninja"
	BUILD_BACKEND="-GNinja"
else
	echo "failed to find ninja"
	BUILD_BACKEND=""
fi

CC=$(which clang)
CXX=$(which clang++)

mkdir "build"
mkdir "build/release"
mkdir "build/debug"

cd "build/release"
CC=$CC CXX=$CXX cmake "-DCMAKE_BUILD_TYPE=Release" "-DDEV_ENVIRONMENT=ON" "-DBUILD_TESTS=ON" "$BUILD_BACKEND"  "../.."
cd "../.."

cd "build/debug"
CC=$CC CXX=$CXX cmake "-DCMAKE_BUILD_TYPE=Debug" "-DDEV_ENVIRONMENT=ON" "-DBUILD_TESTS=ON" "$BUILD_BACKEND"  "../.."
cd "../.."

