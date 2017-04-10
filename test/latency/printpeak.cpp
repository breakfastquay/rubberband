
#include <sndfile.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

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

    size_t maxi = 0, mini = 0, chunkAbsMaxi = 0;
    float max = 0.f, min = 0.f, chunkAbsMax = 0.f;
    int chunkCount = 0;
    int chunkSize = nf/3;

    for (size_t i = 0; i < nf; ++i) {
	float f = ff[i];
	if (f > max) {
	    max = f;
	    maxi = i;
	}
	if (f < min) {
	    min = f;
	    mini = i;
	}
	float af = fabsf(f);
	if (af > 0.05 &&
	    (i == 0 || af > fabsf(ff[i-1])) &&
	    (i == nf-1 || af > fabsf(ff[i+1]))) {
	    cout << "notable peak " << f << " @ " << i << endl;
	}
	if (af > chunkAbsMax) {
	    chunkAbsMax = af;
	    chunkAbsMaxi = i;
	}
	if (i % chunkSize == chunkSize-1) {
	    if (chunkAbsMax > 0.f) {
		cout << "peak in chunk " << chunkCount << " is " << chunkAbsMax
		     << " @ " << chunkAbsMaxi << endl;
	    }
	    chunkCount++;
	    chunkAbsMax = 0;
	    chunkAbsMaxi = i;
	}
    }

    if (chunkAbsMax > 0.f) {
	cout << "peak in chunk " << chunkCount << " is " << chunkAbsMax
	     << " @ " << chunkAbsMaxi << endl;
    }
    cout << "max " << max << " @ " << maxi << endl;
    cout << "min " << min << " @ " << mini << endl;
    sf_close(sndfile);
    
    return 0;
}
