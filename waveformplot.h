#ifndef WAVEFORMPLOT_H
#define WAVEFORMPLOT_H
#include <mgl2/mgl.h>
#include <mgl2/qt.h>
#include <waveform.h>

class WaveFormPlot: public mglDraw
{
    WaveForm m_wave;
    float m_beginSec;
    float m_endSec;
public:
    WaveFormPlot(WaveForm wave);
    int Draw(mglGraph *gr);

    void PlotWaveForm(float begSec=0, float endSec=0);
};

#endif // WAVEFORMPLOT_H
