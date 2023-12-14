
package com.breakfastquay.rubberband.test;

import com.breakfastquay.rubberband.RubberBandStretcher;

import java.util.TreeMap;

public class RubberBandTest
{

    public static void main(String[] args) {

        int channels = 1;
        int rate = 44100;
        
        RubberBandStretcher stretcher = new RubberBandStretcher
            (rate,
             channels,
             RubberBandStretcher.OptionEngineFiner +
             RubberBandStretcher.OptionProcessOffline,
             1.0,
             1.0);

        stretcher.setTimeRatio(1.5);
        stretcher.setPitchScale(0.8);
        
        System.err.println
            (String.format("Channel count: %d\n" +
                           "Time ratio: %f\n" +
                           "Pitch scale: %f\n" +
                           "Preferred start pad: %d\n" +
                           "Start delay: %d\n" +
                           "Process size limit: %d",
                           stretcher.getChannelCount(),
                           stretcher.getTimeRatio(),
                           stretcher.getPitchScale(),
                           stretcher.getPreferredStartPad(),
                           stretcher.getStartDelay(),
                           stretcher.getProcessSizeLimit()
                ));

        int blocksize = 1024;
        int blocks = 400;
        double freq = 440.0;
        
        stretcher.setMaxProcessSize(blocksize);

        TreeMap<Long, Long> keyFrameMap = new TreeMap<Long, Long>();
        keyFrameMap.put((long)(3 * rate), (long)(4 * rate));
        keyFrameMap.put((long)(5 * rate), (long)(5 * rate));
        stretcher.setKeyFrameMap(keyFrameMap);

        float[][] buffer = new float[channels][blocksize];

        int i0 = 0;

        for (int block = 0; block < blocks; ++block) {

            for (int c = 0; c < channels; ++c) {
                for (int i = 0; i < blocksize; ++i) {
                    buffer[c][i] = (float)Math.sin
                        ((double)i0 * freq * Math.PI * 2.0 / (double)rate);
                    if (i0 % rate == 0) {
                        buffer[c][i] = 1.f;
                    }
                    ++i0;
                }
            }
            
            stretcher.study(buffer, block + 1 == blocks);
        }

        i0 = 0;

        for (int block = 0; block < blocks; ++block) {

            for (int c = 0; c < channels; ++c) {
                for (int i = 0; i < blocksize; ++i) {
                    buffer[c][i] = (float)Math.sin
                        ((double)i0 * freq * Math.PI * 2.0 / (double)rate);
                    if (i0 % rate == 0) {
                        buffer[c][i] = 1.f;
                    }
                    ++i0;
                }
            }
            
            stretcher.process(buffer, block + 1 == blocks);

            while (true) {
                int available = stretcher.available();
                if (available <= 0) {
                    break;
                }
                int requested = available;
                if (requested > blocksize) {
                    requested = blocksize;
                }
                int obtained = stretcher.retrieve(buffer, 0, requested);
                for (int i = 0; i < obtained; ++i) {
                    System.out.println(Float.toString(buffer[0][i]));
                }
            }
        }
        
        stretcher.dispose();
    }
    
}

