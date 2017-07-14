#include "autocorrelator.h"

AutoCorrelator::AutoCorrelator(int slidingwindowsize, int analysisFrames)
{
    m_slidingWindowSize=slidingwindowsize;
    m_slidingWindow=new float[slidingwindowsize];
    m_analysisFrames=analysisFrames;
    m_autocorrelationResult=new float[m_analysisFrames];

}

void AutoCorrelator::runAnalysis(WaveForm &input, int beginframe)
{
    float* buf=input.getBuffer();
    for(int i=beginframe;i<input.getBufferSize()-1;i++)
        if(buf[i-1]<0.0&&buf[i+1]>0.0)
        {
            beginframe=i;
            break;
        }
        else if(buf[i-1]>0.0&&buf[i+1]<0.0)
        {
            beginframe=i;
            break;
        }

    for(int i=0;i<m_slidingWindowSize;i++)
        m_slidingWindow[i]=buf[beginframe+i];
    float mean1;
    float mean2;



    float maxval=1e3;
    int maxid=0;
    float f;

    float lastval=1e2;
    float minfreq=1e3;
    float diffsum=0;
    int bestbeg=0;
    int bestend=0;
    for(int i=0;i<m_analysisFrames;i++)
    {
        float m2;

        int beg=i+beginframe;
        int covwinsize=min(m_analysisFrames-i,m_slidingWindowSize);
        int end=beg+covwinsize;
        float cov=covar(m_slidingWindow,0,covwinsize,
                        buf,beg,end,&mean1,&mean2);
        float varData=covar(buf,beg,end,
                            buf,beg,end,&m2,&m2);
        float varWindow=covar(m_slidingWindow,0,covwinsize,m_slidingWindow,0,covwinsize,&mean1,&mean1);
        float cond=(1-cov*cov/(varWindow*varData))*varData;

        float diff=0;

//        mean1=0;
//        mean2=0;
        for(int j=beg;j<end;j++)
            diff+=((m_slidingWindow[j-beg]-mean1)-(buf[j]-mean2)*cov/varData);//*cov/varData;
        diffsum+=diff;
//        diff/=(diffsum*(end-beg));
//        int mid=i+beginframe;
//        int end=m_analysisFrames+beginframe;
//        int beg1=beginframe;
//        int end1=beginframe+i;
//        int end2=beginframe+m_analysisFrames;
//        int beg2=end2-i;
//        float cov=covar(buf,beg1,end1,
//                        buf,beg2,end2,&m1,&m2);
//        float varData=covar(buf,beg2,end2,
//                            buf,beg2,end2);
//        varWindow=covar(buf,beg1,end1,
//                        buf,beg1,end1);

        m_autocorrelationResult[i]=diff*diff/(cond+1e-5);
//        if(i>m_slidingWindowSize)//
            if(m_autocorrelationResult[i]>lastval)
        {
            if(m_autocorrelationResult[i]<maxval//&&fabs(cov/varWindow)<=1.0
                    )//&&input.getSamplerate()/((float)i)<minfreq)
            {
                maxval=m_autocorrelationResult[i];
                maxid=m_analysisFrames-i;
                f=cov/varWindow;
//                if(fabs(f)>1)
//                   f=1.0/f;
                bestbeg=beg;
                bestend=end;
                minfreq=input.getSamplerate()/((float)i);
                cout<<"cond: "<<cond<<endl;
            }
            lastval=-1;
        }
        else
            lastval=m_autocorrelationResult[i];



    }
//    f=0.981011;

    for(int i=maxid+beginframe;i>=beginframe;i--)
    {

        if((buf[i+1]>0&&buf[i-1]<0))
        {
            maxid=i-beginframe;
            break;
        }
        else if((buf[i+1]<0&&buf[i-1]>0))
        {
            maxid=i-beginframe;
            break;
        }
    }
    int overlapwsize=1;
    float covover=covar(buf,beginframe,beginframe+overlapwsize,
                        buf,beginframe+maxid,beginframe+maxid+overlapwsize);
    float varBegin=covar(buf,beginframe,beginframe+overlapwsize,
                         buf,beginframe,beginframe+overlapwsize);
    float varOverlap=covar(buf,beginframe+maxid,beginframe+maxid+overlapwsize,
                         buf,beginframe+maxid,beginframe+maxid+overlapwsize);
    f=0;//covover/(varBegin);
    int count=0;
    for(int i=0;i<overlapwsize;i++)
        if(buf[beginframe+i]!=0)
        {
            f+=(buf[i+beginframe+maxid]-buf[i+beginframe+maxid-1])/(buf[i+beginframe]-buf[i+beginframe-1]);
            count++;
        }
    f/=count;
//    f=1;
    m_baseFrequency=input.getSamplerate()/((float)maxid);
    cout<<"Found base frequency: "<<input.getSamplerate()/((float)maxid)<<endl;
    cout<<"Base framesize: "<<maxid<<endl;
    cout<<"Got transformfactor: "<<f<<endl;
    cout<<"Mean1: "<<mean1<<"Mean2: "<<mean2<<endl;
    int repetitions=1;//abs(log(0.0000001)/log(fabs(f)));
    cout<<"Repetitions: "<<repetitions<<endl;
//    maxid*=4;
    WaveForm repwave(maxid*repetitions,1,input.getSamplerate());
    m_periodFrames=maxid;
    float* repbuf=repwave.getBuffer();
    int pos=0;

//    f=buf[beginframe+maxid]/buf[beginframe];
    for(int r=0;r<repetitions;r++)
        for(int i=0;i<maxid;i++)
        {
            float covover=covar(buf,beginframe,beginframe+overlapwsize,
                                buf,beginframe+i,beginframe+i+overlapwsize);
            float varOverlap=covar(buf,beginframe+i,beginframe+i+overlapwsize,
                                 buf,beginframe+i,beginframe+i+overlapwsize);
//            f=covover/varOverlap;
//            float sd=maxid*repetitions/3.0f;//+(((float)rand())/RAND_MAX)*1000.0;
            float sd =maxid;
            repbuf[pos]=buf[beginframe+i]*pow(f*2,r);//*pow(f/fabs(f),r)*pow(fabs(f),((float)pos/sd));
//            repbuf[pos]=buf[i+beginframe]*exp(-f*pos/sd);
//            if(r)
//                repbuf[pos]+=(mean2)*pow(f,((float)r))-mean1*pow(f,((float)r)-1.0);
            pos++;
        }
    m_repeatedWave=repwave;
    mglQT gr(this,"MathGL examples");
    gr.Run();
}

float AutoCorrelator::mean(const float *buf, int begin, int end)
{
    float res=0;
    for(int i=begin;i<end;i++)
        res+=buf[i];
    res/=(end-begin);
    return res;
}

float AutoCorrelator::covar(const float *buf, int begin, int end, const float *buf2, int begin2, int end2,float* m1,float* m2)
{
    assert(end-begin==end2-begin2);
    float mean1=mean(buf,begin,end);
    float mean2=mean(buf2,begin2,end2);
    int size=end-begin;
    float res=0;
    for(int i=0;i<size;i++)
    {
        res+=(buf[i+begin]-mean1)*(buf2[i+begin2]-mean2);
    }
    res/=(size-1);
    if(m1)
        *m1=mean1;
    if(m2)
        *m2=mean2;
    return res;
}

int AutoCorrelator::Draw(mglGraph *gr)
{

    mglData dat(m_analysisFrames);	// data to for plotting

    float maxval=0;
    float minval=1e3;
    for(long i=0;i<(m_analysisFrames);i++)
    {


        dat.a[i]=m_autocorrelationResult[i];

        maxval=max(maxval,m_autocorrelationResult[i]);
        minval=min(minval,m_autocorrelationResult[i]);


    }

    gr->SetRange('y',minval,maxval);

    gr->Plot(dat);

}
