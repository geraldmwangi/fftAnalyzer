#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <sndfile.hh>
#include <string>
#include <iostream>
#include <assert.h>
using namespace std;
class WaveForm
{
    float* m_buffer;
    int m_buffersize;
    int m_framesize;
    int m_channels;
    int m_samplerate;
    string m_file;
public:
    WaveForm(int frames=0,int channels=0, int samplerate=0);
    WaveForm(string file);
    void writeFile(string fname);
    void readFile();
    float* getBuffer() const
    {
        return m_buffer;
    }

    int getBufferSize() const
    {
        return m_buffersize;
    }
    int getFramesize() const
    {
        return m_framesize;
    }
    int getSamplerate() const
    {
        return m_samplerate;
    }
    int getChannels() const
    {
        return m_channels;
    }

    void resample(float resampleratio);

    void rectifiy();

    WaveForm operator +(const WaveForm& src) const;
    WaveForm &operator=(const WaveForm& src);

};

#endif // WAVEFORM_H
