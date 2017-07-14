#include <iostream>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>
#include <sndfile.h>
#include <waveform.h>
#include <waveformplot.h>
#include <fftanalysis.h>
#include <frequencysieve.h>
#include <autocorrelator.h>
using namespace std;





int main(int argc, char *argv[])
{

//    LogDraw foo;
//    mglQT gr(&foo,"MathGL examples");
//    return gr.Run();

//    WaveForm wave("/home/gerald/Programming/PolyphonicGuitar/Samples/fret4_6_eb_high_norm.wav");
//    wave.readFile();
////    wave.resample(2);
//    WaveForm sum("/home/gerald/Programming/PolyphonicGuitar/Samples/fret4_6_eb_sum_high_norm_eq1k.wav");
//    sum.readFile();
////    wave2.resample(2);
//    WaveForm sum("/home/gerald/Programming/PolyphonicGuitar/Samples/fret4_6_eb_high_norm.wav");
//    WaveForm sum("/home/gerald/Programming/PolyphonicGuitar/Samples/fret4_6_eb_sum_high_norm_eq1k.wav");//=wave+wave2;
//    sum.readFile();
//     WaveForm sum("/home/gerald/Programming/PolyphonicGuitar/Samples/Gerald/fret0_dgbe.wav");
//     sum.readFile();
         WaveForm sum("/home/gerald/Programming/PolyphonicGuitar/Samples/Gerald/3string_fret11.wav");
         sum.readFile();
    WaveForm sum1("/home/gerald/Programming/PolyphonicGuitar/Samples/joeboy/twostring-fret20_eb.wav");
    sum1.readFile();

    WaveForm sum2("/home/gerald/Programming/PolyphonicGuitar/Samples/joeboy/twostring-fret10_eb.wav");
    sum2.readFile();
//    WaveForm sum;
//    sum=sum1+sum2;
    bool useAutoCorr=false;
    int length=512*32;
    float sr=((float)sum.getSamplerate())/(length);
    float begsec=0.2;
    if(useAutoCorr)
    {
        AutoCorrelator corr(1*512,32*512);
        corr.runAnalysis(sum,0.5*sum.getSamplerate());
        sum=corr.m_repeatedWave;
        sum.writeFile("/home/gerald/Programming/PolyphonicGuitar/Samples/joeboy/twostring-fret_chuckrepeat.wav");
        sr=corr.m_baseFrequency;
        length=sum.getBufferSize();
        begsec=0;
    }
    sum.writeFile("/home/gerald/Programming/PolyphonicGuitar/Samples/sumSamples.wav");
    //    sum.rectifiy();
//    sum.resample(1/4.0);
         cout<<"Num samples: "<<sum.getBufferSize()<<endl;
    WaveFormPlot plotwave(sum);

    plotwave.PlotWaveForm(begsec,begsec+((float)length)/sum.getSamplerate());//sum.getBufferSize()/sum.getSamplerate());//(begsec,begsec+512.0f/sum.getSamplerate()*2);

//    FFTAnalysis analysis(sum,corr.m_periodFrames);//sum.getBufferSize());
        FFTAnalysis analysis(sum,length,sr);
    analysis.runAnalysis(begsec);

//    FFTAnalysis analysis1(wave);
//    analysis1.runAnalysis(begsec,4000);

//    FFTAnalysis analysis2(wave2);
//    analysis2.runAnalysis(begsec,4000);

    float thresh=6;
    float potfactor=1.5;
//    FrequencySieve sieve0(analysis1.getPeakFrequencies(),analysis1.getPeakFrequencies(),thresh);
//    cout<<"Computed F0 "<<sieve0.computedF0()<<endl;


//    FrequencySieve sieve1(analysis2.getPeakFrequencies(),analysis2.getPeakFrequencies(),thresh);
//    cout<<"Computed F0 "<<sieve1.computedF0()<<endl;
    srand(time(0));
    FrequencySieve sievesum(analysis.getPeakFrequencies(),analysis.getPeakFrequencies(),thresh);
    cout<<"Sum Computed F0 "<<sievesum.computedF0()<<endl;
    sievesum.printInliers();

    srand(time(0));
    FrequencySieve sievesum2(sievesum.getOutlierFrequencies(),sievesum.getOutlierFrequencies(),thresh*potfactor);
    cout<<"Sum Computed F0 "<<sievesum2.computedF0()<<endl;
    sievesum2.printInliers();
    srand(time(0));
    FrequencySieve sievesum3(sievesum2.getOutlierFrequencies(),sievesum2.getOutlierFrequencies(),thresh*pow(potfactor,2.0));
    cout<<"Sum Computed F0 "<<sievesum3.computedF0()<<endl;
    sievesum3.printInliers();
    srand(time(0));
    FrequencySieve sievesum4(sievesum3.getOutlierFrequencies(),sievesum3.getOutlierFrequencies(),thresh*pow(potfactor,3.0));
    cout<<"Sum Computed F0 "<<sievesum4.computedF0()<<endl;
    sievesum4.printInliers();
}
