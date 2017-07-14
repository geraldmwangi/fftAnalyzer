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
//void FFTAnalysis::medianFilterD(double *array,int arraysize, int filtersize, float power,float* maxval)
//{
//    int peaks[10000];

//    int numpeaks=0;
//    float diff=0;
//    for(long i=1;i<(arraysize-1);i+=1)
//        if(((array[i]-array[i-1])>0.0&&(array[i]-array[i+1])>0.0&&numpeaks<10000)

//                //&&(powerspec[i-1]>(powerspec[i-2]))&&(powerspec[i+1]>(powerspec[i+2]))
//                //&&(powerspec[i-2]>(powerspec[i-3]))&&(powerspec[i+2]>(powerspec[i+3]))
//                )
//        {
////            if(powerspec[i]>0.1&&(m_numpeaks<100))
//            {
//                peaks[numpeaks]=i;
//                diff+=peaks[numpeaks]-peaks[numpeaks-1];
//                numpeaks++;

//            }
//        }
//    diff/=numpeaks;
////    filtersize/=diff;
//    int m=(filtersize-1)/2;
//    float* window=new float[filtersize];
//    double* outarray=new double[arraysize];
//        float timebase=m_frames/((float)(m_wave.getSamplerate()));
//        int delta=1;//timebase*m_baseFrequency;//82.5843;
//    float minval=0;
//    memset(outarray,0,arraysize*sizeof(float));
//    for(int i=0;i<numpeaks;i+=delta)
//    {
//        window[0]=pow(array[peaks[i]],power);
//        for(int j=-m;j<m;j++)
//        {
//            window[j+m]=0;
//            int id=j*delta+i;
//            float x=1;//j*j*2.0/(float)(m*m);
//            if(id>=0&&id<numpeaks)
//                window[j+m]=pow(array[peaks[id]],power);

//        }
//        if(filtersize>1)
//            qsort(window,filtersize,sizeof(float),compare);
//        outarray[peaks[i]]=pow(array[peaks[i]],power)-(window[m]);
//        if(m==0)
//            outarray[peaks[i]]=(window[m]);
////        cout<<pow(array[i],power)<<" "<<window[m]<<endl;



//    }
//    for(int i=0;i<m;i++)
//    outarray[i]=0;
//    double* medianarray=new double[arraysize];
//    memcpy(medianarray,outarray,sizeof(double)*arraysize);
//    qsort(medianarray,arraysize,sizeof(float),compare);
//    minval=0;//medianarray[arraysize/2];
////    for(int i=1;i<arraysize;i++)
////        minval+=(float)outarray[i];//min(minval,(float)outarray[i]);
////    minval/=arraysize;
//        *maxval=-1e3;
//    for(int i=0;i<numpeaks;i++)
//    {
////        outarray[i]=outarray[i]-minval;
//        *maxval=max(*maxval,(float)outarray[peaks[i]]);
//    }
//    for(int i=0;i<numpeaks;i++)
//        outarray[peaks[i]]=(outarray[peaks[i]]-minval)/(*maxval-minval); //(outarray[i]-minval)/(*maxval-minval);
//    *maxval=1;
//    minval=0;
////    for(int i=1;i<arraysize;i++)
////        minval+=(float)outarray[i];//min(minval,(float)outarray[i]);
////    minval/=arraysize;
////    for(int i=0;i<arraysize;i++)
////    {
////        outarray[i]=outarray[i]-minval;
////        *maxval=max(*maxval,(float)outarray[i]);
////    }
////    for(int i=0;i<arraysize;i++)
////        outarray[i]=outarray[i]/(*maxval);
////    *maxval=1;
//    memcpy(array,outarray,sizeof(double)*arraysize);
//    delete [] outarray;
//    delete [] window;
//    delete [] medianarray;
//}
void FFTAnalysis::medianFilterD(double *array,int arraysize, int filtersize, float power,float* maxval)
{
    int m=(filtersize-1)/2;
    float* window=new float[filtersize];
    double* outarray=new double[arraysize];
        float timebase=m_frames/((float)(m_wave.getSamplerate()));
        int delta=1;//timebase*m_baseFrequency;//82.5843;
    float minval=0;
    memset(outarray,0,arraysize*sizeof(float));
    for(int i=0;i<arraysize;i+=delta)
    {
        window[0]=pow(array[i],power);
        for(int j=-m;j<m;j++)
        {
            window[j+m]=0;
            int id=j*delta+i;
            float x=1;//j*j*2.0/(float)(m*m);
            if(id>=0&&id<arraysize)
                window[j+m]=pow(array[id],power);

        }
        if(filtersize>1)
            qsort(window,filtersize,sizeof(float),compare);
        outarray[i]=pow(array[i],power)-(window[m]);
        if(m==0)
            outarray[i]=(window[m]);
//        cout<<pow(array[i],power)<<" "<<window[m]<<endl;



    }
    for(int i=0;i<m;i++)
    outarray[i]=0;
    double* medianarray=new double[arraysize];
    memcpy(medianarray,outarray,sizeof(double)*arraysize);
    qsort(medianarray,arraysize,sizeof(float),compare);
    minval=medianarray[arraysize/2];
//    for(int i=1;i<arraysize;i++)
//        minval+=(float)outarray[i];//min(minval,(float)outarray[i]);
//    minval/=arraysize;
        *maxval=-1e3;
    for(int i=0;i<arraysize;i++)
    {
//        outarray[i]=outarray[i]-minval;
        *maxval=max(*maxval,(float)outarray[i]);
    }
    for(int i=0;i<arraysize;i++)
        outarray[i]=(outarray[i]-minval)/(*maxval-minval); //(outarray[i]-minval)/(*maxval-minval);
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
//void FFTAnalysis::computePowerSpectrumPeaks()
//{
//    int size=m_frames/2.0-1.0;
//    double* powerspec=new double[size];
//    float maxval=0;
//    for(long i=1;i<(size);i++)
//    {

//        float real=m_out[i];
//        float imag=m_out[m_frames-i];
//        float val=((real*real+imag*imag));
//        powerspec[i]=val;
//        maxval=max(val,maxval);


//    }
//    float real=m_out[0];

//    powerspec[0]=(real*real);
//    medianFilterD(powerspec,size,1,0.001f,&maxval);
//    //    medianFilterD(dat.a,size,1,0.001f,&maxval);
////        medianFilterD(powerspec,size,1,1,&maxval);
////    medianFilterD(powerspec,size,5,0.1f,&maxval);
////    medianFilterD(powerspec,size,5,0.1f,&maxval);
//    m_numpeaks=0;
//    float timebase=m_frames/((float)(m_wave.getSamplerate()));
//    int delta=timebase*m_baseFrequency;//82.5843;
//    for(long i=0;i<(size-delta);i+=delta)
//        if((powerspec[i]-(powerspec[i-delta]))>0.01&&(powerspec[i]-(powerspec[i+delta]))>0.01

//                //&&(powerspec[i-1]>(powerspec[i-2]))&&(powerspec[i+1]>(powerspec[i+2]))
//                //&&(powerspec[i-2]>(powerspec[i-3]))&&(powerspec[i+2]>(powerspec[i+3]))
//                )
//        {
//            if(fabs(powerspec[i])>0.2&&(m_numpeaks<100))
//            {
//                m_peaks[m_numpeaks]=i;
//                m_numpeaks++;
//            }
//        }
//    delete [] powerspec;

//}
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
    medianFilterD(powerspec,size,111,0.01f,&maxval);
    //    medianFilterD(dat.a,size,1,0.001f,&maxval);
//        medianFilterD(powerspec,size,1,1,&maxval);
//    medianFilterD(powerspec,size,5,0.1f,&maxval);
//    medianFilterD(powerspec,size,5,0.1f,&maxval);
    m_numpeaks=0;
    int peaks[10000];

    int numpeaks=0;
    for(long i=1;i<(size-1);i+=1)
        if(((powerspec[i]-powerspec[i-1])>0.0&&(powerspec[i]-powerspec[i+1])>0.0&&numpeaks<10000)

                //&&(powerspec[i-1]>(powerspec[i-2]))&&(powerspec[i+1]>(powerspec[i+2]))
                //&&(powerspec[i-2]>(powerspec[i-3]))&&(powerspec[i+2]>(powerspec[i+3]))
                )
        {
//            if(powerspec[i]>0.1&&(m_numpeaks<100))
            {
                peaks[numpeaks]=i;
                numpeaks++;
            }
        }
    for(long i=1;i<(numpeaks-1);i+=1)
        if((powerspec[peaks[i]]>(powerspec[peaks[i-1]]))
                &&(powerspec[peaks[i]]>(powerspec[peaks[i+1]]))

                //&&(powerspec[i-1]>(powerspec[i-2]))&&(powerspec[i+1]>(powerspec[i+2]))
                //&&(powerspec[i-2]>(powerspec[i-3]))&&(powerspec[i+2]>(powerspec[i+3]))
                )
        {
            if(powerspec[peaks[i]]>0.5&&(m_numpeaks<100))
            {
                cout<<powerspec[peaks[i-1]]<<"<"<<powerspec[peaks[i]]
                                          <<">"<<powerspec[peaks[i+1]]<<endl;
                m_peaks[m_numpeaks]=peaks[i];
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

    float windowfreq=((float)m_multiple)/timebase;
    for(long i=1;i<(size);i++)
    {

        float real=m_out[i];
        float imag=m_out[m_frames-i];
//        float r2=pow(real,2.0)-pow(imag,2.0);
//        float im2=2*real*imag;

//        real=r2;
//        imag=im2;
        float val=((real*real+imag*imag));
        float freq=i/timebase;
        float diff=freq-round(freq/windowfreq)*windowfreq;
//        if(diff==0)
//            val=0;
        dat.a[i]=val;
        maxval=max(maxval,val);
        datpeak.a[i]=0;

        xaxis.a[i]=i/timebase;

    }
    float real=m_out[0];

    dat.a[0]=(real*real);
    medianFilterD(dat.a,size,111,0.01,&maxval);
//    medianFilterD(dat.a,size,1,1,&maxval);
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

FFTAnalysis::FFTAnalysis(WaveForm wave, int frames, float baseFrequency)
{
    m_baseFrequency=baseFrequency;
    m_wave=wave;
    m_multiple=1;
    m_out=fftwf_alloc_real(frames*m_multiple);
    m_in=fftwf_alloc_real(frames*m_multiple);
    m_plan=fftwf_plan_r2r_1d(frames*m_multiple,m_in,m_out,FFTW_R2HC, FFTW_ESTIMATE);
    m_frames=frames*m_multiple;
    m_hann_window  = new float[m_frames];
    double sum = 0.0;

    for (uint32_t i = 0; i < m_frames; ++i) {
        float pos=i-frames/2.0;
        float sd=frames/4.0;
//        m_hann_window[i] = exp(-pow(pos/sd,2.0)/2.0);//
        m_hann_window[i]=0.5f - (0.5f * (float) cos (2.0f * M_PI * (float)i / (float)(m_frames)));
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
    int actframes=m_frames/m_multiple;
    int begfr=(m_multiple-1)/2*actframes;
    int endfr=beg+actframes;
    for(int i=0;i<m_frames;i++)
    {
        //if(i>=begfr&&i<endfr)
        {
            int pos=i%actframes;
            float sd=m_frames/300.0f;//+(((float)rand())/RAND_MAX)*100.0;
            m_in[i]=m_hann_window[i]*(buf[pos+beg]);//*exp(-(float)i/sd);//+0.1*((float)rand())/RAND_MAX)/maxval;
        }
//        else
//            m_in[i]=0.0;
//        m_in[i]=buf[pos+beg];
        //cout<<"in: "<<m_in[i]<<" pos:"<<pos<<endl;
    }
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
        if(m_peaks[i]/timebase>100&&m_peaks[i]/timebase<5000)
        {
            res.push_back((m_peaks[i]/timebase));
            cout<<"Get freq: "<<(m_peaks[i]/timebase)<<endl;
        }
    return res;
}
