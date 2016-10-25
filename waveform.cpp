#include "waveform.h"
#include <samplerate.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

WaveForm::WaveForm(int frames, int channels, int samplerate)
{
    m_buffer=0;
    m_buffersize=frames*channels;
    m_samplerate=samplerate;
    m_channels=channels;
    m_framesize=frames;

    if(m_buffersize)
        m_buffer=new float[m_buffersize];
}

WaveForm::WaveForm(string file)
{
    m_file=file;
    m_buffer=0;
    m_buffersize=0;
    m_samplerate=0;
    m_channels=0;
    m_framesize=0;
}

void WaveForm::writeFile(string fname)
{
    SndfileHandle outfile(fname,SFM_WRITE);
//    outfile.
}

void WaveForm::readFile()
{

    SndfileHandle infile( m_file, SFM_READ );
    m_framesize=infile.frames();
    m_channels=infile.channels();

    m_buffersize=infile.frames() * infile.channels();
    m_samplerate=infile.samplerate();
    cout<<"Reading File: "<<m_file<<endl;
    cout<<"Frames: "<<m_framesize<<", Channels: "<<m_channels
       <<" Samplerate: "<<m_samplerate<<endl;
    m_buffer=new float[m_buffersize];
    infile.read( m_buffer , m_buffersize );
    assert(m_channels==1);
}

void WaveForm::resample(float resampleratio)
{
    SRC_DATA data;
    int newframessize=m_framesize*resampleratio;
    int newbuffsize=m_buffersize*resampleratio;
    int newsr=m_samplerate*resampleratio;
    float* newBuffer=new float[newbuffsize];

    data.data_in  = m_buffer;
    data.data_out = newBuffer;

    data.input_frames = m_framesize;
    data.output_frames = newframessize;

    data.end_of_input = 0;
    data.src_ratio = resampleratio;

    int ret = src_simple ( &data, SRC_SINC_BEST_QUALITY, 1 );
    if(ret==0)
    {
        cout<<"Resampling ok"<<endl;
        m_framesize=newframessize;
        m_buffersize=newbuffsize;
        m_samplerate=newsr;
        delete [] m_buffer;
        m_buffer=newBuffer;
    }
    else
    {
        cout<<"Resampling bad!"<<endl;
        delete [] newBuffer;
    }

}

void WaveForm::rectifiy()
{
    for(int i=0;i<m_buffersize;i++)
        m_buffer[i]=fabs(m_buffer[i]);
}

WaveForm WaveForm::operator +(const WaveForm &src) const
{
    assert(m_channels==src.getChannels());
    assert(m_samplerate==src.getSamplerate());
    int newframes=min(m_framesize,src.getFramesize());
    WaveForm sum(newframes,m_channels,m_samplerate);
    float *buf=sum.getBuffer();
    float *srcbuf=src.getBuffer();
    for(int i=0;i<sum.getBufferSize();i++)
        buf[i]=m_buffer[i]+srcbuf[i];
    return sum;
}

WaveForm &WaveForm::operator=(const WaveForm &src)
{
    m_buffersize=src.getBufferSize();
    m_channels=src.getChannels();
    m_framesize=src.getFramesize();
    m_samplerate=src.getSamplerate();
    if(m_buffer)
        delete [] m_buffer;
    m_buffer=new float[m_buffersize];
    memcpy(m_buffer,src.getBuffer(),sizeof(float)*m_buffersize);
    return *this;
}
