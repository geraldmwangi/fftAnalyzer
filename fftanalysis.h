#ifndef FFTANALYSIS_H
#define FFTANALYSIS_H

#include <waveform.h>
#include <fftw3.h>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>
#include <vector>
class FFTAnalysis: public mglDraw
{
    WaveForm m_wave;
    fftwf_plan m_plan;
    float* m_out;
    float *m_in;
    int m_frames;
    int m_peaks[1000];
    int m_numpeaks;
    int m_multiple;
    float* m_hann_window;
    void medianFilter(float* array, int arraysize, int filtersize, float power);
public:
    void computePowerSpectrumPeaks();
    int Draw(mglGraph *gr);
    FFTAnalysis(WaveForm wave,int frames);
    ~FFTAnalysis();
    void runAnalysis(float beginSec);
    std::vector<float> getPeakFrequencies();
private:
    void medianFilterD(double *array, int arraysize, int filtersize, float power,float* maxval);
};

#endif // FFTANALYSIS_H
