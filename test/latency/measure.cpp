
#include <sndfile.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Timing measurement for output of non-transient-preserving mode.
//
// We know that our file contains two impulses, one near the start (or
// at least within the first third of the file) and the other toward
// the end (or at least within the last two-thirds).
//
// These impulses are likely to be smeared, so we want to isolate them
// and find their "middle", i.e. the half-way point between where the
// smeared transient rises out of the noise floor and where it sinks
// back in again. This is not necessarily the same as the peak.
//
// Having located the rough middle, we then look for a peak within the
// area of the middle (the smeared impulse can be asymmetric). This is
// the peak value within the middle 1/4 of the impulse range.

int findTransientCentre(const vector<float> &ff, int i0, int i1)
{
    int i, j, k;
    float threshold = 1e-4f;

    for (i = i0; i < i1; ++i) {
	if (fabsf(ff[i]) > threshold) {
	    cout << "starts at " << i << endl;
	    break;
	}
    }

    if (i == i1) {
	cout << "WARNING: failed to find transient between " << i0
	     << " and " << i1 << endl;
	return i0;
    }
    
    for (j = i1; j >= i0; ) {
	--j;
	if (fabsf(ff[j]) > threshold) {
	    cout << "ends at " << j << endl;
	    break;
	}
    }

    int middle = (i + j) / 2;

    int k0 = i + (j - i) / 4;
    int k1 = j - (j - i) / 4;

    if (k1 <= k0) {
	cout << "WARNING: k1 (" << k1 << ") <= k0 (" << k0 << ")" << endl;
	return k0;
    }

    int peak = i0;
    float peakValue = 0.f;
    
    for (k = k0; k < k1; ++k) {
	if (fabsf(ff[k]) > peakValue) {
	    peakValue = fabsf(ff[k]);
	    peak = k;
	}
    }

    cout << "i0 = " << i0 << ", i1 = " << i1 << ", i = " << i << ", j = " << j << ", middle = " << middle << ", k0 = " << k0 << ", k1 = " << k1 << ", peak = " << peak << ", peakValue = " << peakValue << endl;

    return peak;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
	cerr << "usage: " << argv[0] << " <file.wav>" << endl;
	return 2;
    }

    SF_INFO sfinfo;
    SNDFILE *sndfile = sf_open(argv[1], SFM_READ, &sfinfo);
    if (!sndfile) {
	cerr << "ERROR: Failed to open input file \"" << argv[1] << "\": "
	     << sf_strerror(sndfile) << endl;
	return 1;
    }

    if (sfinfo.channels > 1) {
	cerr << "Mono only please" << endl;
	return 1;
    }

    size_t nf = sfinfo.frames;
    int rate = sfinfo.samplerate;
    
    vector<float> ff(nf, 0.f);
    sf_readf_float(sndfile, &ff[0], nf);

    // before this we seek the first transient, after it we seek the second:
    int division = nf/3;

    int t1 = findTransientCentre(ff, 0, division);
    int t2 = findTransientCentre(ff, division, nf);

    cout << "transient 1 centre @ " << t1 << endl;
    cout << "transient 2 centre @ " << t2 << endl;

    sf_close(sndfile);
    
    return 0;
}
