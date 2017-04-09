#!/bin/bash 

set -eu

( cd ../.. ; make )

if [ ! -f in.wav ]; then
    flac -d in.flac
fi

#sox dirac.wav up.wav pad 100000s 99500s
#sox -v -0.3 dirac.wav down.wav pad 1000s 995s
#sox -m up.wav down.wav testfile.wav

g++ printpeak.cpp -o printpeak -lsndfile
g++ -I../.. printlatency.cpp -o printlatency ../../lib/librubberband.a -lfftw3 -lfftw3f -lsamplerate -lpthread

mkdir -p output

(
    
for timeratio in 0.2 0.4 0.5 0.8 0.999 1.0 1.001 1.2 2.0 2.2 3.4 10.0 ; do
#    for pitchshift in -13 -5 0 5 13 ; do
    for pitchshift in 0 ; do
	#	for rt in N Y ; do
	#	for rt in N; do
	for rt in Y; do
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

		    echo -n "[fftsize $fftsize, in incr $inincr, out incr $outincr] "

		    peak1=$(./printpeak "$outfile" | grep chunk | head -1 | awk '{ print $8; }')
		    peak2=$(./printpeak "$outfile" | grep chunk | tail -1 | awk '{ print $8; }')
		    
		    exp1=$(echo 1000 "$timeratio" '*' p | dc | sed 's/[.].*$//')
		    exp2=$(echo 100000 "$timeratio" '*' p | dc | sed 's/[.].*$//')

		    err1=$(($peak1 - $exp1))
		    err2=$(($peak2 - $exp2))

		    abs1=$(echo "$err1" | sed 's/^-//')
		    abs2=$(echo "$err2" | sed 's/^-//')
		    
		    if [ "$abs1" -lt 3 ]; then
			echo -n "OK ($peak1) "
		    else
			err=$(($peak1 - $exp1))
			echo -n "FAIL (exp $exp1, got $peak1, err $err) "
		    fi

		    if [ "$abs2" -lt 3 ]; then
			echo "OK ($peak2)"
		    else
			err=$(($peak2 - $exp2))
			echo "FAIL (exp $exp2, got $peak2, err $err)"
		    fi

		    rm output/log.txt
		done
	    done
	done
    done
done

) | tee test.log
