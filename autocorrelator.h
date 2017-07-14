#ifndef AUTOCORRELATOR_H
#define AUTOCORRELATOR_H
#include <waveform.h>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>
class AutoCorrelator: public mglDraw
{
    int m_slidingWindowSize;
    float* m_slidingWindow;
    int m_analysisFrames;
    float *m_autocorrelationResult;

public:
    AutoCorrelator(int slidingwindowsize, int analysisFrames);
    void runAnalysis(WaveForm& input, int beginframe);
    float mean(const float* buf,int begin,int end);
    float covar(const float* buf, int begin, int end,
                const float* buf2, int begin2, int end2, float *m1=0, float *m2=0);
    int Draw(mglGraph *gr);
    WaveForm m_repeatedWave;
    int m_periodFrames;
    float m_baseFrequency;
};

#endif // AUTOCORRELATOR_H
