#!/bin/bash
set -eu
if [ ! -d /Applications ]; then
    make -f otherbuilds/Makefile.linux
    g++ main/main.cpp lib/librubberband.a -I. -Isrc -o test -lsndfile -lsamplerate -lpthread
    ./test -V
    g++ main/main.cpp single/RubberBandSingle.cpp -o test_single -lsndfile
    ./test_single -V
else
    make -f otherbuilds/Makefile.macos
    c++ main/main.cpp lib/librubberband.a -I. -Isrc -o test -lsndfile -lsamplerate -framework Accelerate
    ./test -V
    make -f otherbuilds/Makefile.macos clean
    make -f otherbuilds/Makefile.ios
    c++ main/main.cpp single/RubberBandSingle.cpp -o test_single -lsndfile
    ./test_single -V
fi
