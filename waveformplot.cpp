#include "waveformplot.h"

WaveFormPlot::WaveFormPlot(WaveForm wave)
{
    m_wave=wave;
}

int WaveFormPlot::Draw(mglGraph *gr)
{
    float* buf=m_wave.getBuffer();

    int sr=m_wave.getSamplerate();
    int end=(m_endSec>0)?sr*m_endSec:m_wave.getBufferSize();
    end=min(end,m_wave.getBufferSize());
    int beg=m_beginSec*sr;
    mglData dat(end-beg);	// data to for plotting
    mglData xdat(end-beg);
    for(long i=0;i<(end-beg);i++)
    {
        int pos=i+beg;
        dat.a[i]=buf[pos];
        xdat.a[i]=((float)pos)/sr;
    }
    //dat.Set(buf,m_wave.getBufferSize());
    //        mglGraph gr;		// class for plot drawing

    gr->SetRange('x',xdat);
    gr->Plot(dat);		// plot surface

    gr->Axis();			// draw axi
}

void WaveFormPlot::PlotWaveForm(float begSec,float endSec)
{
    m_beginSec=begSec;
    m_endSec=endSec;
    mglQT gr(this,"MathGL examples");
    gr.Run();
}
