#include "fftanalysis.h"
#include <xtract/libxtract.h>
#include <cstdlib>
int compare( const void* a, const void* b)
{
     float int_a = * ( (float*) a );
     float int_b = * ( (float*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}
void FFTAnalysis::medianFilter(float *array,int arraysize, int filtersize, float power)
{
    int m=(filtersize-1)/2;
    float* window=new float[filtersize];
    float* outarray=new float[arraysize];
    float maxval=0;
    float minval=1;
    for(int i=0;i<arraysize;i++)
    {
        window[0]=pow(array[i],power);
        for(int j=-m;j<m;j++)
        {
            window[j+m]=0;
            int id=j+i;
            if(id>=0&&id<arraysize)
                window[j+m]=pow(array[id],power);

        }
        qsort(window,filtersize,sizeof(float),compare);
        outarray[i]=(window[m]);
        maxval=max(maxval,(float)outarray[i]);


    }
    for(int i=1;i<arraysize;i++)
        minval=min(minval,(float)outarray[i]);
    for(int i=0;i<arraysize;i++)
        outarray[i]=(outarray[i]-minval)/(maxval-minval);
    maxval=1;
    memcpy(array,outarray,sizeof(float)*arraysize);
    delete [] outarray;
    delete [] window;
}
void FFTAnalysis::medianFilterD(double *array,int arraysize, int filtersize, float power,float* maxval)
{
    int m=(filtersize-1)/2;
    float* window=new float[filtersize];
    double* outarray=new double[arraysize];
    *maxval=0;
    float minval=0;
    for(int i=0;i<arraysize;i++)
    {
        window[0]=pow(array[i],power);
        for(int j=-m;j<m;j++)
        {
            window[j+m]=0;
            int id=j+i;
            float x=1;//j*j*2.0/(float)(m*m);
            if(id>=0&&id<arraysize)
                window[j+m]=pow(array[id]*exp(-x),power);

        }
        if(filtersize>1)
            qsort(window,filtersize,sizeof(float),compare);
        outarray[i]=(window[m]);



    }
    double* medianarray=new double[arraysize];
    memcpy(medianarray,outarray,sizeof(double)*arraysize);
    qsort(medianarray,arraysize,sizeof(float),compare);
    minval=medianarray[arraysize/2];
//    for(int i=1;i<arraysize;i++)
//        minval+=(float)outarray[i];//min(minval,(float)outarray[i]);
//    minval/=arraysize;
    for(int i=0;i<arraysize;i++)
    {
        outarray[i]=outarray[i]-minval;
        *maxval=max(*maxval,(float)outarray[i]);
    }
    for(int i=0;i<arraysize;i++)
        outarray[i]=outarray[i]/(*maxval); //(outarray[i]-minval)/(*maxval-minval);
    *maxval=1;
    minval=0;
//    for(int i=1;i<arraysize;i++)
//        minval+=(float)outarray[i];//min(minval,(float)outarray[i]);
//    minval/=arraysize;
//    for(int i=0;i<arraysize;i++)
//    {
//        outarray[i]=outarray[i]-minval;
//        *maxval=max(*maxval,(float)outarray[i]);
//    }
//    for(int i=0;i<arraysize;i++)
//        outarray[i]=outarray[i]/(*maxval);
//    *maxval=1;
    memcpy(array,outarray,sizeof(double)*arraysize);
    delete [] outarray;
    delete [] window;
    delete [] medianarray;
}
void FFTAnalysis::computePowerSpectrumPeaks()
{
    int size=m_frames/2.0-1.0;
    double* powerspec=new double[size];
    float maxval=0;
    for(long i=1;i<(size);i++)
    {

        float real=m_out[i];
        float imag=m_out[m_frames-i];
        float val=((real*real+imag*imag));
        powerspec[i]=val;
        maxval=max(val,maxval);


    }
    float real=m_out[0];

    powerspec[0]=(real*real);
    medianFilterD(powerspec,size,1,0.001f,&maxval);
//    medianFilterD(powerspec,size,5,0.1f,&maxval);
//    medianFilterD(powerspec,size,5,0.1f,&maxval);
    m_numpeaks=0;
    for(long i=3;i<(size-1);i+=1)
        if((powerspec[i]>(powerspec[i-1]))&&(powerspec[i]>(powerspec[i+1]))

                //&&(powerspec[i-1]>(powerspec[i-2]))&&(powerspec[i+1]>(powerspec[i+2]))
                //&&(powerspec[i-2]>(powerspec[i-3]))&&(powerspec[i+2]>(powerspec[i+3]))
                )
        {
            if(powerspec[i]>0.5&&(m_numpeaks<100))
            {
                m_peaks[m_numpeaks]=i;
                m_numpeaks++;
            }
        }
    delete [] powerspec;

}

int FFTAnalysis::Draw(mglGraph *gr)
{
    float timebase=m_frames/((float)(m_wave.getSamplerate()));
    int size=4000*timebase;//m_frames/2.0-1.0;
    mglData dat(size);	// data to for plotting
    mglData datpeak(size);
    mglData xaxis(size);
    mglData xpeak(m_numpeaks);
    mglData ypeak(m_numpeaks);
    float maxval=0;

    for(long i=1;i<(size);i++)
    {

        float real=m_out[i];
        float imag=m_out[m_frames-i];
//        float r2=pow(real,2.0)-pow(imag,2.0);
//        float im2=2*real*imag;

//        real=r2;
//        imag=im2;
        float val=((real*real+imag*imag));
        dat.a[i]=val;
        maxval=max(maxval,val);
        datpeak.a[i]=0;

        xaxis.a[i]=i/timebase;

    }
    float real=m_out[0];

    dat.a[0]=(real*real);
    medianFilterD(dat.a,size,1,0.001f,&maxval);
//    medianFilterD(dat.a,size,5,0.1f,&maxval);
//    medianFilterD(dat.a,size,5,0.1f,&maxval);
//    medianFilterD(dat.a,size,5,0.1,&maxval);
//    medianFilterD(dat.a,size,3,0.1,&maxval);

    //dat.Set(buf,m_wave.getBufferSize());
    //        mglGraph gr;		// class for plot drawing


            // plot surface

   // gr->Box();
    gr->SetRange('y',0.0,maxval);
    gr->SetRange('x',0.0,size/timebase);
    gr->Plot(xaxis,dat);
    gr->Axis();
    for(int i=0;i<m_numpeaks&&m_peaks[i]<size;i++)
    {
        datpeak.a[m_peaks[i]]=maxval;

    }
//    gr->Plot(datpeak,".");
    //gr->Plot(ypeak,xpeak,"+");
  //
  //  gr->Axis(); gr->Grid("xy","g;");
   // gr->Axis();			// draw axi
}

FFTAnalysis::FFTAnalysis(WaveForm wave, int frames)
{
    m_wave=wave;
    m_out=fftwf_alloc_real(frames);
    m_in=fftwf_alloc_real(frames);
    m_plan=fftwf_plan_r2r_1d(frames,m_in,m_out,FFTW_R2HC, FFTW_ESTIMATE);
    m_frames=frames;
    m_hann_window  = new float[m_frames];
    double sum = 0.0;

    for (uint32_t i = 0; i < m_frames; ++i) {
        m_hann_window[i] = 0.5f - (0.5f * (float) cos (2.0f * M_PI * (float)i / (float)(m_frames)));
        sum += m_hann_window[i];
    }
    const double isum = 2.0 / sum;
    for (uint32_t i = 0; i < m_frames; ++i) {
        m_hann_window[i] *= isum;
    }
}

FFTAnalysis::~FFTAnalysis()
{
    if(m_out)
        fftwf_free(m_out);
    if(m_in)
        fftwf_free(m_in);
    fftwf_destroy_plan(m_plan);
    delete [] m_hann_window;
}

void FFTAnalysis::runAnalysis(float beginSec)
{
    float timebase=m_frames/((float)(m_wave.getSamplerate()));
    cout<<"Freq res:"<<1.0f/timebase<<endl;
    int beg=m_wave.getSamplerate()*beginSec;
    int end=beg+m_frames;






    float* buf=m_wave.getBuffer();
    float maxval=1;
//    maxval=0;
//    for(int i=0;i<m_frames;i++)
//        maxval=max(buf[i+beg],maxval);
    for(int i=0;i<m_frames;i++)
        m_in[i]=m_hann_window[i]*buf[i+beg]/maxval;
    fftwf_execute(m_plan);
    computePowerSpectrumPeaks();
    mglQT gr(this,"MathGL examples");
    gr.Run();


}

std::vector<float> FFTAnalysis::getPeakFrequencies()
{
    std::vector<float> res;
    computePowerSpectrumPeaks();
    float timebase=m_frames/((float)(m_wave.getSamplerate()));
    for(int i=0;i<m_numpeaks;i++)
        if(m_peaks[i]/timebase>500)
        {
            res.push_back(m_peaks[i]/timebase);
           // cout<<"Get freq: "<<m_peaks[i]/timebase<<endl;
        }
    return res;
}
