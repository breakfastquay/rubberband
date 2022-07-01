#include "finer/R3StretcherImpl.h"

int main(int argc, char **argv)
{
    RubberBand::R3StretcherImpl::Parameters parameters(44100.0, 2);
    RubberBand::R3StretcherImpl impl(parameters);
}

