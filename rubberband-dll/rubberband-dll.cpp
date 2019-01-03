// rubberband-dll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"

using namespace RubberBand;

#define EXPORT extern "C" __declspec(dllexport)

#if _WIN64
#define CONVENTION __stdcall
#else
#define CONVENTION __cdecl
#endif

#define O(rbs) reinterpret_cast<RubberBandStretcher *>(rbs)

EXPORT HANDLE CONVENTION RubberBandStretcher_Create(
	size_t sampleRate,
	size_t channels,
	int options = RubberBandStretcher::DefaultOptions,
	double initialTimeRatio = 1.0,
	double initialPitchScale = 1.0)
{
	return (HANDLE)new RubberBandStretcher(sampleRate, channels, options, initialTimeRatio, initialPitchScale);
}

EXPORT void CONVENTION RubberBandStretcher_Delete(HANDLE *rbs)
{
	delete O(rbs);
}

EXPORT void CONVENTION RubberBandStretcher_Reset(HANDLE *rbs)
{
	O(rbs)->reset();
}

EXPORT void CONVENTION RubberBandStretcher_SetTimeRatio(HANDLE *rbs, double ratio)
{
	O(rbs)->setTimeRatio(ratio);
}

EXPORT void CONVENTION RubberBandStretcher_SetPitchScale(HANDLE *rbs, double scale)
{
	O(rbs)->setPitchScale(scale);
}

EXPORT double CONVENTION RubberBandStretcher_GetTimeRatio(HANDLE *rbs)
{
	return O(rbs)->getTimeRatio();
}

EXPORT double CONVENTION RubberBandStretcher_GetPitchScale(HANDLE *rbs)
{
	return O(rbs)->getPitchScale();
}

EXPORT size_t CONVENTION RubberBandStretcher_GetLatency(HANDLE *rbs)
{
	return O(rbs)->getLatency();
}

EXPORT void CONVENTION RubberBandStretcher_SetTransientsOption(HANDLE *rbs, int options)
{
	O(rbs)->setTransientsOption(options);
}

EXPORT void CONVENTION RubberBandStretcher_SetDetectorOption(HANDLE *rbs, int options)
{
	O(rbs)->setDetectorOption(options);
}

EXPORT void CONVENTION RubberBandStretcher_SetPhaseOption(HANDLE *rbs, int options)
{
	O(rbs)->setPhaseOption(options);
}

EXPORT void CONVENTION RubberBandStretcher_SetFormantOption(HANDLE *rbs, int options)
{
	O(rbs)->setFormantOption(options);
}

EXPORT void CONVENTION RubberBandStretcher_SetPitchOption(HANDLE *rbs, int options)
{
	O(rbs)->setPitchOption(options);
}

EXPORT void CONVENTION RubberBandStretcher_SetExpectedInputDuration(HANDLE *rbs, size_t samples)
{
	O(rbs)->setExpectedInputDuration(samples);
}

EXPORT void CONVENTION RubberBandStretcher_SetMaxProcessSize(HANDLE *rbs, size_t samples)
{
	O(rbs)->setMaxProcessSize(samples);
}

EXPORT size_t CONVENTION RubberBandStretcher_GetSamplesRequired(HANDLE *rbs)
{
	return O(rbs)->getSamplesRequired();
}

EXPORT void CONVENTION RubberBandStretcher_SetKeyFrameMap(HANDLE *rbs, size_t *mappingData, int numberOfMappings)
{
	std::map<size_t, size_t> map;

	for (int i = 0; i < numberOfMappings; i++)
		map[mappingData[i + i]] = mappingData[i + i + 1];

	O(rbs)->setKeyFrameMap(map);
}

EXPORT void CONVENTION RubberBandStretcher_Study(HANDLE *rbs, float *input_flat, size_t samples, int channels, bool final)
{
	std::vector<float *> input;

	for (int i = 0; i < channels; i++)
		input.push_back(&input_flat[samples * i]);

	O(rbs)->study(&input[0], samples, final);
}

EXPORT void CONVENTION RubberBandStretcher_Process(HANDLE *rbs, float *input_flat, size_t samples, int channels, bool final)
{
	std::vector<float *> input;

	for (int i = 0; i < channels; i++)
		input.push_back(&input_flat[samples * i]);

	O(rbs)->process(&input[0], samples, final);
}

EXPORT int CONVENTION RubberBandStretcher_Available(HANDLE *rbs)
{
	return O(rbs)->available();
}

EXPORT size_t CONVENTION RubberBandStretcher_Retrieve(HANDLE *rbs, float *output_flat, size_t samples, int channels)
{
	std::vector<float *> output;

	for (int i = 0; i < channels; i++)
		output.push_back(&output_flat[samples * i]);

	return O(rbs)->retrieve(&output[0], samples);
}

EXPORT float CONVENTION RubberBandStretcher_GetFrequencyCutoff(HANDLE *rbs, int n)
{
	return O(rbs)->getFrequencyCutoff(n);
}

EXPORT void CONVENTION RubberBandStretcher_SetFrequencyCutoff(HANDLE *rbs, int n, float f)
{
	O(rbs)->setFrequencyCutoff(n, f);
}

EXPORT size_t CONVENTION RubberBandStretcher_GetInputIncrement(HANDLE *rbs)
{
	return O(rbs)->getInputIncrement();
}

EXPORT size_t CONVENTION RubberBandStretcher_GetOutputIncrements(HANDLE *rbs, int *buffer, size_t bufferSize)
{
	std::vector<int> increments = O(rbs)->getOutputIncrements();

	for (size_t i = 0; (i < bufferSize) && (i < increments.size()); i++)
		buffer[i] = increments[i];

	return increments.size();
}

EXPORT size_t CONVENTION RubberBandStretcher_GetPhaseResetCurve(HANDLE *rbs, float *buffer, size_t bufferSize)
{
	std::vector<float> curve = O(rbs)->getPhaseResetCurve();

	for (size_t i = 0; (i < bufferSize) && (i < curve.size()); i++)
		buffer[i] = curve[i];

	return curve.size();
}

EXPORT size_t CONVENTION RubberBandStretcher_GetExactTimePoints(HANDLE *rbs, int *buffer, size_t bufferSize)
{
	std::vector<int> timePoints = O(rbs)->getExactTimePoints();

	for (size_t i = 0; (i < bufferSize) && (i < timePoints.size()); i++)
		buffer[i] = timePoints[i];

	return timePoints.size();
}

EXPORT size_t CONVENTION RubberBandStretcher_GetChannelCount(HANDLE *rbs)
{
	return O(rbs)->getChannelCount();
}

EXPORT void CONVENTION RubberBandStretcher_CalculateStretch(HANDLE *rbs)
{
	O(rbs)->calculateStretch();
}

EXPORT void CONVENTION RubberBandStretcher_SetDebugLevel(HANDLE *rbs, int level)
{
	O(rbs)->setDebugLevel(level);
}

EXPORT void CONVENTION RubberBandStretcher_SetDefaultDebugLevel(int level)
{
	RubberBandStretcher::setDefaultDebugLevel(level);
}
