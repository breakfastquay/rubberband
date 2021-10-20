/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2021 Particular Programs Ltd.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.

    Alternatively, if you have a valid commercial licence for the
    Rubber Band Library obtained by agreement with the copyright
    holders, you may redistribute and/or modify it under the terms
    described in that licence.

    If you wish to distribute code using the Rubber Band Library
    under terms other than those of the GNU General Public License,
    you must obtain a valid commercial licence before doing so.
*/

#include "../rubberband/RubberBandStretcher.h"

#include <iostream>
#include <sndfile.h>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <cstring>
#include <string>

#include <fstream>

#include "../src/system/sysutils.h"

#ifdef _MSC_VER
#include "../src/getopt/getopt.h"
#else
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#endif

#include "../src/base/Profiler.h"

#ifdef _WIN32
using RubberBand::gettimeofday;
#endif

#ifdef _MSC_VER
using RubberBand::usleep;
#define strdup _strdup
#endif

using RubberBand::RubberBandStretcher;

using std::cerr;
using std::endl;

double tempo_convert(const char *str)
{
    char *d = strchr((char *)str, ':');

    if (!d || !*d) {
        double m = atof(str);
        if (m != 0.0) return 1.0 / m;
        else return 1.0;
    }

    char *a = strdup(str);
    char *b = strdup(d+1);
    a[d-str] = '\0';
    double m = atof(a);
    double n = atof(b);
    free(a);
    free(b);
    if (n != 0.0 && m != 0.0) return m / n;
    else return 1.0;
}

int main(int argc, char **argv)
{
    int c;

    double ratio = 1.0;
    double duration = 0.0;
    double pitchshift = 0.0;
    double frequencyshift = 1.0;
    int debug = 0;
    bool realtime = false;
    bool precise = true;
    int threading = 0;
    bool lamination = true;
    bool longwin = false;
    bool shortwin = false;
    bool smoothing = false;
    bool hqpitch = false;
    bool formant = false;
    bool together = false;
    bool crispchanged = false;
    int crispness = -1;
    bool help = false;
    bool version = false;
    bool quiet = false;

    bool haveRatio = false;

    std::string timeMapFile;
    std::string freqMapFile;
    std::string pitchMapFile;
    bool freqOrPitchMapSpecified = false;

    enum {
        NoTransients,
        BandLimitedTransients,
        Transients
    } transients = Transients;

    enum {
        CompoundDetector,
        PercussiveDetector,
        SoftDetector
    } detector = CompoundDetector;

    bool ignoreClipping = false;

    while (1) {
        int optionIndex = 0;

        static struct option longOpts[] = {
            { "help",          0, 0, 'h' },
            { "version",       0, 0, 'V' },
            { "time",          1, 0, 't' },
            { "tempo",         1, 0, 'T' },
            { "duration",      1, 0, 'D' },
            { "pitch",         1, 0, 'p' },
            { "frequency",     1, 0, 'f' },
            { "crisp",         1, 0, 'c' },
            { "crispness",     1, 0, 'c' },
            { "debug",         1, 0, 'd' },
            { "realtime",      0, 0, 'R' },
            { "loose",         0, 0, 'L' },
            { "precise",       0, 0, 'P' },
            { "formant",       0, 0, 'F' },
            { "no-threads",    0, 0, '0' },
            { "no-transients", 0, 0, '1' },
            { "no-lamination", 0, 0, '2' },
            { "centre-focus",  0, 0, '7' },
            { "window-long",   0, 0, '3' },
            { "window-short",  0, 0, '4' },
            { "bl-transients", 0, 0, '8' },
            { "detector-perc", 0, 0, '5' },
            { "detector-soft", 0, 0, '6' },
            { "smoothing",     0, 0, '9' },
            { "pitch-hq",      0, 0, '%' },
            { "threads",       0, 0, '@' },
            { "quiet",         0, 0, 'q' },
            { "timemap",       1, 0, 'M' },
            { "freqmap",       1, 0, 'Q' },
            { "pitchmap",      1, 0, 'C' },
            { "ignore-clipping", 0, 0, 'i' },
            { 0, 0, 0, 0 }
        };

        c = getopt_long(argc, argv,
                        "t:p:d:RLPFc:f:T:D:qhVM:",
                        longOpts, &optionIndex);
        if (c == -1) break;

        switch (c) {
        case 'h': help = true; break;
        case 'V': version = true; break;
        case 't': ratio *= atof(optarg); haveRatio = true; break;
        case 'T': ratio *= tempo_convert(optarg); haveRatio = true; break;
        case 'D': duration = atof(optarg); haveRatio = true; break;
        case 'p': pitchshift = atof(optarg); haveRatio = true; break;
        case 'f': frequencyshift = atof(optarg); haveRatio = true; break;
        case 'd': debug = atoi(optarg); break;
        case 'R': realtime = true; break;
        case 'L': precise = false; break;
        case 'P': precise = true; break;
        case 'F': formant = true; break;
        case '0': threading = 1; break;
        case '@': threading = 2; break;
        case '1': transients = NoTransients; crispchanged = true; break;
        case '2': lamination = false; crispchanged = true; break;
        case '3': longwin = true; crispchanged = true; break;
        case '4': shortwin = true; crispchanged = true; break;
        case '5': detector = PercussiveDetector; crispchanged = true; break;
        case '6': detector = SoftDetector; crispchanged = true; break;
        case '7': together = true; break;
        case '8': transients = BandLimitedTransients; crispchanged = true; break;
        case '9': smoothing = true; crispchanged = true; break;
        case '%': hqpitch = true; break;
        case 'c': crispness = atoi(optarg); break;
        case 'q': quiet = true; break;
        case 'M': timeMapFile = optarg; break;
        case 'Q': freqMapFile = optarg; freqOrPitchMapSpecified = true; break;
        case 'C': pitchMapFile = optarg; freqOrPitchMapSpecified = true; break;
        case 'i': ignoreClipping = true; break;
        default:  help = true; break;
        }
    }

    if (version) {
        cerr << RUBBERBAND_VERSION << endl;
        return 0;
    }

    if (freqOrPitchMapSpecified) {
        if (freqMapFile != "" && pitchMapFile != "") {
            cerr << "ERROR: Please specify either pitch map or frequency map, not both" << endl;
            return 1;
        }
        haveRatio = true;
        realtime = true;
    }
    
    if (help || !haveRatio || optind + 2 != argc) {
        cerr << endl;
	cerr << "Rubber Band" << endl;
        cerr << "An audio time-stretching and pitch-shifting library and utility program." << endl;
	cerr << "Copyright 2007-2021 Particular Programs Ltd." << endl;
        cerr << endl;
	cerr << "   Usage: " << argv[0] << " [options] <infile.wav> <outfile.wav>" << endl;
        cerr << endl;
        cerr << "You must specify at least one of the following time and pitch ratio options." << endl;
        cerr << endl;
        cerr << "  -t<X>, --time <X>       Stretch to X times original duration, or" << endl;
        cerr << "  -T<X>, --tempo <X>      Change tempo by multiple X (same as --time 1/X), or" << endl;
        cerr << "  -T<X>, --tempo <X>:<Y>  Change tempo from X to Y (same as --time X/Y), or" << endl;
        cerr << "  -D<X>, --duration <X>   Stretch or squash to make output file X seconds long" << endl;
        cerr << endl;
        cerr << "  -p<X>, --pitch <X>      Raise pitch by X semitones, or" << endl;
        cerr << "  -f<X>, --frequency <X>  Change frequency by multiple X" << endl;
        cerr << endl;
        cerr << "The following options provide ways of making the time and frequency ratios" << endl;
        cerr << "change during the audio." << endl;
        cerr << endl;
        cerr << "  -M<F>, --timemap <F>    Use file F as the source for time map" << endl;
        cerr << endl;
        cerr << "  A time map (or key-frame map) file contains a series of lines, each with two" << endl;
        cerr << "  sample frame numbers separated by a single space. These are source and" << endl;
        cerr << "  target frames for fixed time points within the audio data, defining a varying" << endl;
        cerr << "  stretch factor through the audio. When supplying a time map you must specify" << endl;
        cerr << "  an overall stretch factor using -t, -T, or -D as well, to determine the" << endl;
        cerr << "  total output duration." << endl;
        cerr << endl;
        cerr << "         --pitchmap <F>   Use file F as the source for pitch map" << endl;
        cerr << endl;
        cerr << "  A pitch map file contains a series of lines, each with two values: the input" << endl;
        cerr << "  sample frame number and a pitch offset in semitones, separated by a single" << endl;
        cerr << "  space. These specify a varying pitch factor through the audio. The offsets" << endl;
        cerr << "  are all relative to an initial offset specified by the pitch or frequency" << endl;
        cerr << "  option, or relative to no shift if neither was specified. Offsets are" << endl;
        cerr << "  not cumulative. This option implies realtime mode (-R) and also enables a" << endl;
        cerr << "  high-consistency pitch shifting mode, appropriate for dynamic pitch changes." << endl;
        cerr << "  Because of the use of realtime mode, the overall duration will not be exact." << endl;
        cerr << endl;
        cerr << "         --freqmap <F>    Use file F as the source for frequency map" << endl;
        cerr << endl;
        cerr << "  A frequency map file is like a pitch map, except that its second column" << endl;
        cerr << "  lists frequency multipliers rather than pitch offsets (like the difference" << endl;
        cerr << "  between pitch and frequency options above)." << endl;
        cerr << endl;
        cerr << "The following options provide a simple way to adjust the sound. See below" << endl;
        cerr << "for more details." << endl;
        cerr << endl;
        cerr << "  -c<N>, --crisp <N>      Crispness (N = 0,1,2,3,4,5,6); default 5 (see below)" << endl;
        cerr << "  -F,    --formant        Enable formant preservation when pitch shifting" << endl;
        cerr << endl;
        cerr << "The remaining options fine-tune the processing mode and stretch algorithm." << endl;
        cerr << "These are mostly included for test purposes; the default settings and standard" << endl;
        cerr << "crispness parameter are intended to provide the best sounding set of options" << endl;
        cerr << "for most situations. The default is to use none of these options." << endl;
        cerr << endl;
        cerr << "  -L,    --loose          Relax timing in hope of better transient preservation" << endl;
        cerr << "  -P,    --precise        Ignored: The opposite of -L, this is default from 1.6" << endl;
        cerr << "  -R,    --realtime       Select realtime mode (implies --no-threads)" << endl;
        cerr << "         --no-threads     No extra threads regardless of CPU and channel count" << endl;
        cerr << "         --threads        Assume multi-CPU even if only one CPU is identified" << endl;
        cerr << "         --no-transients  Disable phase resynchronisation at transients" << endl;
        cerr << "         --bl-transients  Band-limit phase resync to extreme frequencies" << endl;
        cerr << "         --no-lamination  Disable phase lamination" << endl;
        cerr << "         --window-long    Use longer processing window (actual size may vary)" << endl;
        cerr << "         --window-short   Use shorter processing window" << endl;
        cerr << "         --smoothing      Apply window presum and time-domain smoothing" << endl;
        cerr << "         --detector-perc  Use percussive transient detector (as in pre-1.5)" << endl;
        cerr << "         --detector-soft  Use soft transient detector" << endl;
        cerr << "         --pitch-hq       In RT mode, use a slower, higher quality pitch shift" << endl;
        cerr << "         --centre-focus   Preserve focus of centre material in stereo" << endl;
        cerr << "                          (at a cost in width and individual channel quality)" << endl;
        cerr << "         --ignore-clipping Ignore clipping at output; the default is to restart" << endl;
        cerr << "                          with reduced gain if clipping occurs" << endl;
        cerr << endl;
        cerr << "  -d<N>, --debug <N>      Select debug level (N = 0,1,2,3); default 0, full 3" << endl;
        cerr << "                          (N.B. debug level 3 includes audible ticks in output)" << endl;
        cerr << "  -q,    --quiet          Suppress progress output" << endl;
        cerr << endl;
        cerr << "  -V,    --version        Show version number and exit" << endl;
        cerr << "  -h,    --help           Show this help" << endl;
        cerr << endl;
        cerr << "\"Crispness\" levels:" << endl;
        cerr << "  -c 0   equivalent to --no-transients --no-lamination --window-long" << endl;
        cerr << "  -c 1   equivalent to --detector-soft --no-lamination --window-long (for piano)" << endl;
        cerr << "  -c 2   equivalent to --no-transients --no-lamination" << endl;
        cerr << "  -c 3   equivalent to --no-transients" << endl;
        cerr << "  -c 4   equivalent to --bl-transients" << endl;
        cerr << "  -c 5   default processing options" << endl;
        cerr << "  -c 6   equivalent to --no-lamination --window-short (may be good for drums)" << endl;
        cerr << endl;
        return 2;
    }

    if (ratio <= 0.0) {
        cerr << "ERROR: Invalid time ratio " << ratio << endl;
        return 1;
    }

    if (crispness >= 0 && crispchanged) {
        cerr << "WARNING: Both crispness option and transients, lamination or window options" << endl;
        cerr << "         provided -- crispness will override these other options" << endl;
    }

    if (hqpitch && freqOrPitchMapSpecified) {
        cerr << "WARNING: High-quality pitch mode selected, but frequency or pitch map file is" << endl;
        cerr << "         provided -- pitch mode will be overridden by high-consistency mode" << endl;
        hqpitch = false;
    }

    switch (crispness) {
    case -1: crispness = 5; break;
    case 0: detector = CompoundDetector; transients = NoTransients; lamination = false; longwin = true; shortwin = false; break;
    case 1: detector = SoftDetector; transients = Transients; lamination = false; longwin = true; shortwin = false; break;
    case 2: detector = CompoundDetector; transients = NoTransients; lamination = false; longwin = false; shortwin = false; break;
    case 3: detector = CompoundDetector; transients = NoTransients; lamination = true; longwin = false; shortwin = false; break;
    case 4: detector = CompoundDetector; transients = BandLimitedTransients; lamination = true; longwin = false; shortwin = false; break;
    case 5: detector = CompoundDetector; transients = Transients; lamination = true; longwin = false; shortwin = false; break;
    case 6: detector = CompoundDetector; transients = Transients; lamination = false; longwin = false; shortwin = true; break;
    };

    if (!quiet) {
        cerr << "Using crispness level: " << crispness << " (";
        switch (crispness) {
        case 0: cerr << "Mushy"; break;
        case 1: cerr << "Piano"; break;
        case 2: cerr << "Smooth"; break;
        case 3: cerr << "Balanced multitimbral mixture"; break;
        case 4: cerr << "Unpitched percussion with stable notes"; break;
        case 5: cerr << "Crisp monophonic instrumental"; break;
        case 6: cerr << "Unpitched solo percussion"; break;
        }
        cerr << ")" << endl;
    }

    std::map<size_t, size_t> timeMap;
    if (timeMapFile != "") {
        std::ifstream ifile(timeMapFile.c_str());
        if (!ifile.is_open()) {
            cerr << "ERROR: Failed to open time map file \""
                 << timeMapFile << "\"" << endl;
            return 1;
        }
        std::string line;
        int lineno = 0;
        while (!ifile.eof()) {
            std::getline(ifile, line);
            while (line.length() > 0 && line[0] == ' ') {
                line = line.substr(1);
            }
            if (line == "") {
                ++lineno;
                continue;
            }
            std::string::size_type i = line.find_first_of(" ");
            if (i == std::string::npos) {
                cerr << "ERROR: Time map file \"" << timeMapFile
                     << "\" is malformed at line " << lineno << endl;
                return 1;
            }
            size_t source = atoi(line.substr(0, i).c_str());
            while (i < line.length() && line[i] == ' ') ++i;
            size_t target = atoi(line.substr(i).c_str());
            timeMap[source] = target;
            if (debug > 0) {
                cerr << "adding mapping from " << source << " to " << target << endl;
            }
            ++lineno;
        }
        ifile.close();

        if (!quiet) {
            cerr << "Read " << timeMap.size() << " line(s) from time map file" << endl;
        }
    }

    std::map<size_t, double> freqMap;

    if (freqOrPitchMapSpecified) {
        std::string file = freqMapFile;
        bool convertFromPitch = false;
        if (pitchMapFile != "") {
            file = pitchMapFile;
            convertFromPitch = true;
        }
        std::ifstream ifile(file.c_str());
        if (!ifile.is_open()) {
            cerr << "ERROR: Failed to open map file \"" << file << "\"" << endl;
            return 1;
        }
        std::string line;
        int lineno = 0;
        while (!ifile.eof()) {
            std::getline(ifile, line);
            while (line.length() > 0 && line[0] == ' ') {
                line = line.substr(1);
            }
            if (line == "") {
                ++lineno;
                continue;
            }
            std::string::size_type i = line.find_first_of(" ");
            if (i == std::string::npos) {
                cerr << "ERROR: Map file \"" << file
                     << "\" is malformed at line " << lineno << endl;
                return 1;
            }
            size_t source = atoi(line.substr(0, i).c_str());
            while (i < line.length() && line[i] == ' ') ++i;
            double freq = atof(line.substr(i).c_str());
            if (convertFromPitch) {
                freq = pow(2.0, freq / 12.0);
            }
            freqMap[source] = freq;
            if (debug > 0) {
                cerr << "adding mapping for source frame " << source << " of frequency multiplier " << freq << endl;
            }
            ++lineno;
        }
        ifile.close();

        if (!quiet) {
            cerr << "Read " << freqMap.size() << " line(s) from frequency map file" << endl;
        }
    }

    char *fileName = strdup(argv[optind++]);
    char *fileNameOut = strdup(argv[optind++]);

    SNDFILE *sndfile;
    SNDFILE *sndfileOut;
    SF_INFO sfinfo;
    SF_INFO sfinfoOut;
    memset(&sfinfo, 0, sizeof(SF_INFO));
    memset(&sfinfoOut, 0, sizeof(SF_INFO));

    sndfile = sf_open(fileName, SFM_READ, &sfinfo);
    if (!sndfile) {
        cerr << "ERROR: Failed to open input file \"" << fileName << "\": "
             << sf_strerror(sndfile) << endl;
        return 1;
    }

    if (sfinfo.samplerate == 0) {
        cerr << "ERROR: File lacks sample rate in header" << endl;
        return 1;
    }

    if (duration != 0.0) {
        if (sfinfo.frames == 0) {
            cerr << "ERROR: File lacks frame count in header, cannot use --duration" << endl;
            return 1;
        }
        double induration = double(sfinfo.frames) / double(sfinfo.samplerate);
        if (induration != 0.0) ratio = duration / induration;
    }
    
    sfinfoOut.channels = sfinfo.channels;
    sfinfoOut.format = sfinfo.format;
    sfinfoOut.frames = int(sfinfo.frames * ratio + 0.1);
    sfinfoOut.samplerate = sfinfo.samplerate;
    sfinfoOut.sections = sfinfo.sections;
    sfinfoOut.seekable = sfinfo.seekable;
    
    sndfileOut = sf_open(fileNameOut, SFM_WRITE, &sfinfoOut) ;
    if (!sndfileOut) {
        cerr << "ERROR: Failed to open output file \"" << fileNameOut << "\" for writing: "
             << sf_strerror(sndfileOut) << endl;
        return 1;
    }

    RubberBandStretcher::Options options = 0;
    if (realtime)    options |= RubberBandStretcher::OptionProcessRealTime;
    if (precise)     options |= RubberBandStretcher::OptionStretchPrecise;
    if (!lamination) options |= RubberBandStretcher::OptionPhaseIndependent;
    if (longwin)     options |= RubberBandStretcher::OptionWindowLong;
    if (shortwin)    options |= RubberBandStretcher::OptionWindowShort;
    if (smoothing)   options |= RubberBandStretcher::OptionSmoothingOn;
    if (formant)     options |= RubberBandStretcher::OptionFormantPreserved;
    if (together)    options |= RubberBandStretcher::OptionChannelsTogether;

    if (freqOrPitchMapSpecified) {
        options |= RubberBandStretcher::OptionPitchHighConsistency;
    } else if (hqpitch) {
        options |= RubberBandStretcher::OptionPitchHighQuality;
    }
    
    switch (threading) {
    case 0:
        options |= RubberBandStretcher::OptionThreadingAuto;
        break;
    case 1:
        options |= RubberBandStretcher::OptionThreadingNever;
        break;
    case 2:
        options |= RubberBandStretcher::OptionThreadingAlways;
        break;
    }

    switch (transients) {
    case NoTransients:
        options |= RubberBandStretcher::OptionTransientsSmooth;
        break;
    case BandLimitedTransients:
        options |= RubberBandStretcher::OptionTransientsMixed;
        break;
    case Transients:
        options |= RubberBandStretcher::OptionTransientsCrisp;
        break;
    }

    switch (detector) {
    case CompoundDetector:
        options |= RubberBandStretcher::OptionDetectorCompound;
        break;
    case PercussiveDetector:
        options |= RubberBandStretcher::OptionDetectorPercussive;
        break;
    case SoftDetector:
        options |= RubberBandStretcher::OptionDetectorSoft;
        break;
    }

    if (pitchshift != 0.0) {
        frequencyshift *= pow(2.0, pitchshift / 12.0);
    }

    cerr << "Using time ratio " << ratio;

    if (!freqOrPitchMapSpecified) {
        cerr << " and frequency ratio " << frequencyshift << endl;
    } else {
        cerr << " and initial frequency ratio " << frequencyshift << endl;
    }
    
#ifdef _WIN32
    RubberBand::
#endif
    timeval tv;
    (void)gettimeofday(&tv, 0);
    
    RubberBandStretcher::setDefaultDebugLevel(debug);

    size_t countIn = 0, countOut = 0;

    float gain = 1.f;
    bool successful = false;
    
    const size_t channels = sfinfo.channels;
    const int bs = 1024;
    
    float **cbuf = new float *[channels];
    for (size_t c = 0; c < channels; ++c) {
        cbuf[c] = new float[bs];
    }
    float *ibuf = new float[channels * bs];

    int thisBlockSize;

    while (!successful) { // we may have to repeat with a modified
                          // gain, if clipping occurs
        successful = true;

        RubberBandStretcher ts(sfinfo.samplerate, channels, options,
                               ratio, frequencyshift);
        ts.setExpectedInputDuration(sfinfo.frames);

        int frame = 0;
        int percent = 0;

        sf_seek(sndfile, 0, SEEK_SET);

        if (!realtime) {

            if (!quiet) {
                cerr << "Pass 1: Studying..." << endl;
            }

            while (frame < sfinfo.frames) {

                int count = -1;
                if ((count = sf_readf_float(sndfile, ibuf, bs)) <= 0) break;
        
                for (size_t c = 0; c < channels; ++c) {
                    for (int i = 0; i < count; ++i) {
                        cbuf[c][i] = ibuf[i * channels + c];
                    }
                }

                bool final = (frame + bs >= sfinfo.frames);

                ts.study(cbuf, count, final);

                int p = int((double(frame) * 100.0) / sfinfo.frames);
                if (p > percent || frame == 0) {
                    percent = p;
                    if (!quiet) {
                        cerr << "\r" << percent << "% ";
                    }
                }

                frame += bs;
            }

            if (!quiet) {
                cerr << "\rCalculating profile..." << endl;
            }

            sf_seek(sndfile, 0, SEEK_SET);
        }

        frame = 0;
        percent = 0;

        if (!timeMap.empty()) {
            ts.setKeyFrameMap(timeMap);
        }

        std::map<size_t, double>::const_iterator freqMapItr = freqMap.begin();
    
        countIn = 0;
        countOut = 0;
        bool clipping = false;

        // The stretcher only pads the start in offline mode; to avoid
        // a fade in at the start, we pad it manually in RT mode
        int toDrop = 0;
        if (realtime) {
            toDrop = int(ts.getLatency());
            int toPad = int(round(toDrop * frequencyshift));
            if (debug > 0) {
                cerr << "padding start with " << toPad
                     << " samples in RT mode, will drop " << toDrop
                     << " at output" << endl;
            }
            if (toPad > 0) {
                for (size_t c = 0; c < channels; ++c) {
                    for (int i = 0; i < bs; ++i) {
                        cbuf[c][i] = 0.f;
                    }
                }
                while (toPad > 0) {
                    int p = toPad;
                    if (p > bs) p = bs;
                    ts.process(cbuf, p, false);
                    toPad -= p;
                }
            }
        }                
        
        while (frame < sfinfo.frames) {

            thisBlockSize = bs;

            while (freqMapItr != freqMap.end()) {
                size_t nextFreqFrame = freqMapItr->first;
                if (nextFreqFrame <= countIn) {
                    double s = frequencyshift * freqMapItr->second;
                    if (debug > 0) {
                        cerr << "at frame " << countIn
                             << " (requested at " << freqMapItr->first
                             << " [NOT] plus latency " << ts.getLatency()
                             << ") updating frequency ratio to " << s << endl;
                    }
                    ts.setPitchScale(s);
                    ++freqMapItr;
                } else {
                    if (nextFreqFrame < countIn + thisBlockSize) {
                        thisBlockSize = nextFreqFrame - countIn;
                    }
                    break;
                }
            }

            int count = -1;
            if ((count = sf_readf_float(sndfile, ibuf, thisBlockSize)) < 0) {
                break;
            }
        
            countIn += count;

            for (size_t c = 0; c < channels; ++c) {
                for (int i = 0; i < count; ++i) {
                    cbuf[c][i] = ibuf[i * channels + c];
                }
            }

            bool final = (frame + thisBlockSize >= sfinfo.frames);

            if (debug > 2) {
                cerr << "count = " << count << ", bs = " << thisBlockSize << ", frame = " << frame << ", frames = " << sfinfo.frames << ", final = " << final << endl;
            }

            ts.process(cbuf, count, final);

            int avail;
            while ((avail = ts.available()) > 0) {
                if (debug > 1) {
                    cerr << "available = " << avail << endl;
                }

                thisBlockSize = avail;
                if (thisBlockSize > bs) {
                    thisBlockSize = bs;
                }
                
                if (toDrop > 0) {
                    int dropHere = toDrop;
                    if (dropHere > thisBlockSize) {
                        dropHere = thisBlockSize;
                    }
                    if (debug > 1) {
                        cerr << "toDrop = " << toDrop << ", dropping "
                             << dropHere << " of " << avail << endl;
                    }
                    ts.retrieve(cbuf, dropHere);
                    toDrop -= dropHere;
                    avail -= dropHere;
                    continue;
                }
                
                if (debug > 2) {
                    cerr << "retrieving block of " << thisBlockSize << endl;
                }
                ts.retrieve(cbuf, thisBlockSize);
                
                if (realtime && final) {
                    // (in offline mode the stretcher handles this itself)
                    size_t ideal = size_t(countIn * ratio);
                    if (debug > 2) {
                        cerr << "at end, ideal = " << ideal
                             << ", countOut = " << countOut
                             << ", thisBlockSize = " << thisBlockSize << endl;
                    }
                    if (countOut + thisBlockSize > ideal) {
                        thisBlockSize = ideal - countOut;
                        if (debug > 1) {
                            cerr << "truncated final block to " << thisBlockSize
                                 << endl;
                        }
                    }
                }
                
                countOut += thisBlockSize;
                
                for (size_t c = 0; c < channels; ++c) {
                    for (int i = 0; i < thisBlockSize; ++i) {
                        float value = gain * cbuf[c][i];
                        if (ignoreClipping) { // i.e. just clamp, don't bail out
                            if (value > 1.f) value = 1.f;
                            if (value < -1.f) value = -1.f;
                        } else {
                            if (value >= 1.f || value < -1.f) {
                                clipping = true;
                                gain = (0.999f / fabsf(cbuf[c][i]));
                            }
                        }
                        ibuf[i * channels + c] = value;
                    }
                }
                sf_writef_float(sndfileOut, ibuf, thisBlockSize);
            }

            if (clipping) {
                if (!quiet) {
                    cerr << "NOTE: Clipping detected at output sample "
                         << countOut << ", restarting with "
                         << "reduced gain of " << gain
                         << " (supply --ignore-clipping to avoid this)" << endl;
                }
                const float mingain = 0.75f;
                if (gain < mingain) {
                    cerr << "WARNING: Clipped values were implausibly high: "
                         << "something wrong with input or process - "
                         << "not reducing gain below " << mingain << endl;
                    gain = mingain;
                    ignoreClipping = true;
                }
                successful = false;
                break;
            }
            
            if (frame == 0 && !realtime && !quiet) {
                cerr << "Pass 2: Processing..." << endl;
            }

            int p = int((double(frame) * 100.0) / sfinfo.frames);
            if (p > percent || frame == 0) {
                percent = p;
                if (!quiet) {
                    cerr << "\r" << percent << "% ";
                }
            }

            frame += count;
        }

        if (!successful) {
            sf_seek(sndfile, 0, SEEK_SET);
            sf_seek(sndfileOut, 0, SEEK_SET);
            continue;
        }
    
        if (!quiet) {
            cerr << "\r    " << endl;
        }

        int avail;
        while ((avail = ts.available()) >= 0) {
            if (debug > 1) {
                cerr << "(completing) available = " << avail << endl;
            }

            if (avail == 0) {
                if (realtime ||
                    (options & RubberBandStretcher::OptionThreadingNever)) {
                    break;
                } else {
                    usleep(10000);
                }
            }
            
            thisBlockSize = avail;
            if (thisBlockSize > bs) {
                thisBlockSize = bs;
            }
                
            ts.retrieve(cbuf, thisBlockSize);

            countOut += thisBlockSize;
                
            for (size_t c = 0; c < channels; ++c) {
                for (int i = 0; i < thisBlockSize; ++i) {
                    float value = gain * cbuf[c][i];
                    if (value > 1.f) value = 1.f;
                    if (value < -1.f) value = -1.f;
                    ibuf[i * channels + c] = value;
                }
            }
                
            sf_writef_float(sndfileOut, ibuf, thisBlockSize);
        }
    }

    delete[] ibuf;

    for (size_t c = 0; c < channels; ++c) {
        delete[] cbuf[c];
    }
    delete[] cbuf;

    sf_close(sndfile);
    sf_close(sndfileOut);

    free(fileName);
    free(fileNameOut);

    if (!quiet) {

        cerr << "in: " << countIn << ", out: " << countOut << ", ratio: " << float(countOut)/float(countIn) << ", ideal output: " << lrint(countIn * ratio) << ", error: " << abs(lrint(countIn * ratio) - int(countOut)) << endl;

#ifdef _WIN32
        RubberBand::
#endif
            timeval etv;
        (void)gettimeofday(&etv, 0);
        
        etv.tv_sec -= tv.tv_sec;
        if (etv.tv_usec < tv.tv_usec) {
            etv.tv_usec += 1000000;
            etv.tv_sec -= 1;
        }
        etv.tv_usec -= tv.tv_usec;
        
        double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
        cerr << "elapsed time: " << sec << " sec, in frames/sec: " << countIn/sec << ", out frames/sec: " << countOut/sec << endl;
    }

    RubberBand::Profiler::dump();
    
    return 0;
}


