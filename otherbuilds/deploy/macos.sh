#!/bin/bash
set -eu
if [ ! -f ../rba/deploy/macos/notarize.sh ]; then
    echo "need notarize script in ../rba/deploy/macos"
fi
version=$(grep '^ *version:' meson.build | head -1 | sed "s/^.*'\([0-9][0-9.]*\)'.*$/\1/")
echo
echo "Packaging command-line utility for Mac for Rubber Band v$version..."
echo
if [ -f /usr/local/lib/libsndfile.dylib ]; then
    echo "(WARNING: libsndfile dynamic library found in /usr/local/lib, be sure that you aren't about to combine this external dependency with the hardened runtime)"
fi
rm -rf build
PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/ meson build --cross-file ./cross/macos-universal.txt
ninja -C build
echo
echo "Check the following version number: it should read $version"
./build/rubberband -V
echo
key="Developer ID Application: Particular Programs Ltd (73F996B92S)"
mkdir -p packages
( cd build
  codesign -s "$key" -fv --options runtime rubberband
  zipfile="rubberband-$version-gpl-executable-macos.zip"
  rm -f "$zipfile"
  ditto -c -k rubberband "$zipfile"
  ../../rba/deploy/macos/notarize.sh "$zipfile" com.breakfastquay.rubberband
)
package_dir="rubberband-$version-gpl-executable-macos"
rm -rf "$package_dir"
mkdir "$package_dir"
cp build/rubberband "$package_dir"
cp CHANGELOG README.md COPYING "$package_dir"
tar cvjf "$package_dir.tar.bz2" "$package_dir"
mv "$package_dir.tar.bz2" packages/
rm -rf "$package_dir"
echo
echo "Done, package is in packages/$package_dir.tar.bz2"

