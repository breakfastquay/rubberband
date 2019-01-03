using System;
using System.Runtime.InteropServices;

namespace RubberBand
{
	internal class RubberBandNativeMethodsx64
	{
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_Create(IntPtr sampleRate, IntPtr channels, int options = 0, double initialTimeRatio = 1.0, double initialPitchScale = 1.0);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_Delete(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_Reset(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetTimeRatio(IntPtr rbs, double ratio);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetPitchScale(IntPtr rbs, double scale);
		[DllImport("rubberband-dll-x64")]
		public static extern double RubberBandStretcher_GetTimeRatio(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern double RubberBandStretcher_GetPitchScale(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetLatency(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetTransientsOption(IntPtr rbs, int options);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetDetectorOption(IntPtr rbs, int options);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetPhaseOption(IntPtr rbs, int options);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetFormantOption(IntPtr rbs, int options);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetPitchOption(IntPtr rbs, int options);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetExpectedInputDuration(IntPtr rbs, IntPtr samples);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetMaxProcessSize(IntPtr rbs, IntPtr samples);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetSamplesRequired(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetKeyFrameMap(IntPtr rbs, IntPtr[] mappingData, int numberOfMappings);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_Study(IntPtr rbs, float[] input, IntPtr samples, int channels, bool final);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_Process(IntPtr rbs, float[] input, IntPtr samples, int channels, bool final);
		[DllImport("rubberband-dll-x64")]
		public static extern int RubberBandStretcher_Available(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_Retrieve(IntPtr rbs, float[] output, IntPtr samples, int channels);
		[DllImport("rubberband-dll-x64")]
		public static extern float RubberBandStretcher_GetFrequencyCutoff(IntPtr rbs, int n);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetFrequencyCutoff(IntPtr rbs, int n, float f);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetInputIncrement(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetOutputIncrements(IntPtr rbs, int[] buffer, IntPtr bufferSize);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetPhaseResetCurve(IntPtr rbs, float[] buffer, IntPtr bufferSize);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetExactTimePoints(IntPtr rbs, int[] buffer, IntPtr bufferSize);
		[DllImport("rubberband-dll-x64")]
		public static extern IntPtr RubberBandStretcher_GetChannelCount(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_CalculateStretch(IntPtr rbs);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetDebugLevel(IntPtr rbs, int level);
		[DllImport("rubberband-dll-x64")]
		public static extern void RubberBandStretcher_SetDefaultDebugLevel(int level);
	}
}
