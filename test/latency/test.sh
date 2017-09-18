#!/bin/bash 

set -eu

( cd ../.. ; make )

rm -f in.wav

#if [ ! -f in.wav ]; then
#    flac -d in.flac
#fi

#sox dirac.wav up.wav pad 100000s 99500s
#sox -v -1.0 dirac.wav down.wav pad 1000s 995s
#sox -m up.wav down.wav testfile.wav
#cp testfile.wav in.wav

sox dirac.wav 1.wav pad 1000s 
sox -v -1.0 dirac.wav 2.wav pad 50000s 
sox dirac.wav 3.wav pad 100000s 50000s
sox -m 1.wav 2.wav 3.wav in.wav

g++ printpeak.cpp -o printpeak -lsndfile
g++ measure.cpp -o measure -lsndfile
g++ -I../.. printlatency.cpp -o printlatency ../../lib/librubberband.a -lfftw3 -lfftw3f -lsamplerate -lpthread

mkdir -p output

(
    
#    for pitchshift in -13 -5 0 5 13 ; do
    for pitchshift in 0 ; do
	#	for rt in N Y ; do
	#	for rt in N; do
#	for rt in Y N; do
	for rt in Y; do
#	    for window in L M S ; do
	    for window in M ; do
#		for pitchhq in N Y ; do
		for pitchhq in N ; do
		    for timeratio in 0.2 0.4 0.5 0.8 0.999 1.0 1.001 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.95 2.0 2.05 2.1 2.2 2.3 2.4 2.5 2.55 2.6 2.7 2.8 2.9 3.0 3.4 4.0 10.0 ; do
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
		    ../../bin/rubberband $rtopt $winopt $pitchhqopt --no-delay-comp --no-transients --no-lamination \
					 --time "$timeratio" \
					 --pitch "$pitchshift" \
					 in.wav \
					 "$outfile" -d1 > output/log.txt 2>&1
		    ../../bin/rubberband $rtopt $winopt $pitchhqopt --no-transients \
					 --time "$timeratio" \
					 --pitch "$pitchshift" \
					 drums.wav \
					 "$outdrums" >/dev/null 2>&1
		    fftsize=$(grep 'fft size =' output/log.txt | head -1 | sed 's/^.*fft size = \([0-9]*\).*$/\1/')
		    inincr=$(grep ', increment =' output/log.txt | head -1 | sed 's/^.*, increment = \([0-9]*\).*$/\1/')
		    outincr=$(grep 'output increment =' output/log.txt | head -1 | sed 's/^.*output increment = \([0-9]*\).*$/\1/')
		    delay=$(grep 'reported output delay =' output/log.txt | head -1 | sed 's/^.*output delay = \([0-9]*\).*$/\1/')

		    echo -n "[fftsize $fftsize, in incr $inincr, out incr $outincr, out delay $delay] "

#		    peak1=$(./printpeak "$outfile" | grep chunk | head -1 | awk '{ print $8; }')
#		    peak2=$(./printpeak "$outfile" | grep chunk | tail -n +2 | head -1 | awk '{ print $8; }')

		    peak1=$(./measure "$outfile" | grep 'transient 1' | awk '{ print $5; }')
		    peak2=$(./measure "$outfile" | grep 'transient 2' | awk '{ print $5; }')
		    peak3=$(./measure "$outfile" | grep 'transient 3' | awk '{ print $5; }')
		    
		    exp1=$(echo 1000 "$timeratio" '*' p | dc | sed 's/[.].*$//')
		    exp2=$(echo 50000 "$timeratio" '*' p | dc | sed 's/[.].*$//')
		    exp3=$(echo 100000 "$timeratio" '*' p | dc | sed 's/[.].*$//')

		    exp1=$(($exp1 + $delay))
		    exp2=$(($exp2 + $delay))
		    exp3=$(($exp3 + $delay))
		    
		    err1=$(($peak1 - $exp1))
		    err2=$(($peak2 - $exp2))
		    err3=$(($peak3 - $exp3))

		    abs1=$(echo "$err1" | sed 's/^-//')
		    abs2=$(echo "$err2" | sed 's/^-//')
		    abs3=$(echo "$err3" | sed 's/^-//')
		    
		    if [ "$abs1" -lt 3 ]; then
			echo -n "OK ($peak1) "
		    else
			err=$(($peak1 - $exp1))
			echo -n "FAIL (exp $exp1, got $peak1, err $err) "
		    fi

		    if [ "$abs2" -lt 3 ]; then
			echo -n "OK ($peak2) "
		    else
			err=$(($peak2 - $exp2))
			echo -n "FAIL (exp $exp2, got $peak2, err $err) "
		    fi

		    if [ "$abs3" -lt 3 ]; then
			echo "OK ($peak3)"
		    else
			err=$(($peak3 - $exp3))
			echo "FAIL (exp $exp3, got $peak3, err $err)"
		    fi

		    rm output/log.txt
		    done
		    echo
	    done
	done
    done
done
    
) | tee test.log
