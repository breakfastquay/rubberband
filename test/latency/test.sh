#!/bin/bash 

set -eu

( cd ../.. ; make )

sox dirac.wav up.wav pad 100000s 99500s
# this doesn't work, whichever signal is louder/earlier ends up having
# both the max and min values at more substantial stretch/squash
# factors
sox -v -0.3 dirac.wav down.wav pad 1000s 995s
sox -m up.wav down.wav testfile.wav

g++ printpeak.cpp -o printpeak -lsndfile
g++ -I../.. printlatency.cpp -o printlatency ../../lib/librubberband.a -lfftw3 -lfftw3f -lsamplerate -lpthread

mkdir -p output

(
    
for timeratio in 0.2 0.4 0.9 1.0 1.001 1.2 2.2 3.4 ; do
#    for pitchshift in -13 -5 0 5 13 ; do
    for pitchshift in 0 ; do
	#	for rt in N Y ; do
	for rt in N; do
#	    for window in L M S ; do
	    for window in M ; do
#		for pitchhq in N Y ; do
		for pitchhq in N ; do
		    echo -n "time $timeratio pitch $pitchshift rt $rt win $window hq $pitchhq -> "
		    rtopt=""
		    case $rt in
			Y) rtopt="--realtime";;
		    esac
		    winopt=""
		    case $window in
			L) winopt="--window-long";;
			S) winopt="--window-short";;
		    esac
		    pitchhqopt=""
		    case $pitchhq in
			Y) pitchhqopt="--pitch-hq";;
		    esac
		    outfile="output/t_${timeratio}_${pitchshift}_R=${rt}_W=${window}_P=${pitchhq}.wav"
		    outdrums="output/d_${timeratio}_${pitchshift}_R=${rt}_W=${window}_P=${pitchhq}.wav"
		    ../../bin/rubberband $rtopt $winopt $pitchhqopt --no-transients \
					 --time "$timeratio" \
					 --pitch "$pitchshift" \
					 testfile.wav \
					 "$outfile" -d1 > output/log.txt 2>&1
		    ../../bin/rubberband $rtopt $winopt $pitchhqopt --no-transients \
					 --time "$timeratio" \
					 --pitch "$pitchshift" \
					 drums.wav \
					 "$outdrums" >/dev/null 2>&1
		    fftsize=$(grep 'fft size =' output/log.txt | head -1 | sed 's/^.*fft size = \([0-9]*\).*$/\1/')
		    inincr=$(grep ', increment =' output/log.txt | head -1 | sed 's/^.*, increment = \([0-9]*\).*$/\1/')
		    outincr=$(grep 'output increment =' output/log.txt | head -1 | sed 's/^.*output increment = \([0-9]*\).*$/\1/')
		    peakpos=$(./printpeak "$outfile" | grep max | awk '{ print $4; }')
		    expected=$(echo 100000 "$timeratio" '*' p | dc | sed 's/[.].*$//')
		    alternate=$(($expected + 1))
		    if [ "$peakpos" = "$expected" ] || [ "$peakpos" = "$alternate" ]; then
			echo "OK ($peakpos)"
		    else
			err=$(($expected - $peakpos))
			err=$(echo "$err" | sed 's/^-//')
			echo "FAIL (exp $expected, got $peakpos, err $err, fftsize $fftsize, in incr $inincr, out incr $outincr)"
		    fi

		    rm output/log.txt
		done
	    done
	done
    done
done

) | tee test.log
