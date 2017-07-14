#include "frequencysieve.h"



FrequencySieve::FrequencySieve(vector<float> freqs, vector<float> amplitudes, double thresh):
    Ransac(1,0.999,thresh)
{
    m_freqs=freqs;
    m_amplitudes=amplitudes;
    if(freqs.size()-1)
    {
        for(int i=0;i<(freqs.size()-1);i++)
            for(int j=(i+1);j<freqs.size();j++)
                m_idxs.push_back(pair<int,int>(i,j));
        addData(m_idxs);
        computeRansac();
    }
}

vector<double> FrequencySieve::errorfunction(pair<int, int> dataid, const vector<float> &currF0, const vector<int> &inliers)
{
    //Are the two frequencies of the duplet multiples of the given f0
    float first=m_freqs[dataid.first];
    float second=m_freqs[dataid.second];
    float firstamp=m_amplitudes[dataid.first];
    float secondamp=m_amplitudes[dataid.second];

    int mult1=round(first/currF0[0]);
    int mult2=round(second/currF0[0]);
    double err=fabs(first-mult1*currF0[0])+fabs(second-mult2*currF0[0]);
    int mult3=1;//round(fabs(first-second)/currF0[0]);
    double err2=fabs(fabs(first-second)-mult3*currF0[0]);//+(secondamp>firstamp)?(secondamp-firstamp):0;
//            err/=3.0;
//            err+=1.0/fabs(first-second);
//    double err=fabs(fabs(first-second)-currF0[0]);//+(secondamp>firstamp)?(secondamp-firstamp):0;
//    if(0)
    double err3=0;
    bool firstin=false;
    bool secin=false;
    for(int i=0;i<inliers.size();i++)
    {
        pair<int,int> inpair=m_data[inliers[i]];

        float fdata1=m_freqs[dataid.first];
        float fdata2=m_freqs[dataid.second];
        float fin1=m_freqs[inpair.first];
        float fin2=m_freqs[inpair.second];

        float eps=m_inlierthresh;///2.0;
        if(fabs(fdata1-fin1)<eps&&dataid.first!=inpair.first)
            err3+=m_inlierthresh*10000;
        if(fabs(fdata2-fin1)<eps&&dataid.second!=inpair.first)
            err3+=m_inlierthresh*10000;
        if(fabs(fdata1-fin2)<eps&&dataid.first!=inpair.second)
            err3+=m_inlierthresh*10000;
        if(fabs(fdata2-fin2)<eps&&dataid.second!=inpair.second)
            err3+=m_inlierthresh*10000;
        if(first==fin1||first==fin2)
            firstin=true;
        if(second==fin1||second==fin2)
            secin=true;



    }

    vector<double> res;
    mult3=fabs(mult1-mult2);//round(fabs(first-second)/currF0[0]);
//    if(mult3>1)
//        res.push_back((fabs(fabs(first-second)-mult3*currF0[0]) +
//                fabs(first-mult1*currF0[0])+fabs(second-mult2*currF0[0]))/3.0);
//    else
//        res.push_back(m_inlierthresh);
    if(firstin&&secin)
        res.push_back(0);
    else
        res.push_back((err+err2)+err3);
//    if(fabs(first-second)==110.25)
//        res.push_back(m_inlierthresh*10);
//    res.push_back(err);
//    res.push_back(err2);
//    res.push_back(err3);
    return res;///norm;

}

vector<float> FrequencySieve::computeParameters(const vector<pair<int, int> > &data, int numunknowns)
{
    vector<float> newF0;
    //We get a list of duplets and compute the f0 as the median of the difference data[i].first-data[i].second
    vector<float>freqs;
    for(int i=0;i<data.size();i++)
        freqs.push_back(fabs(m_freqs[data[i].first]-m_freqs[data[i].second]));//freq+=fabs(data[i].first-data[i].second);//min(fabs(data[i].first-data[i].second),freq);
    sort(freqs.begin(),freqs.end());
    double freq=freqs[0];//[freqs.size()/2];
    if(freqs.size()>1)
    {
        if(freqs.size()%2==0)
            freq=freqs[(freqs.size())/2.0];
        else
            freq=freqs[(freqs.size()+1)/2.0];
    }
    newF0.push_back(freq);

    return newF0;
}

void FrequencySieve::printInliers()
{
    vector<bool> mask(m_freqs.size(),false);
    vector<int > inl=inliers();
    for(int i=0;i<inl.size();i++)
    {
        pair<int,int> in=m_idxs[inl[i]];
        mask[in.first]=true;
        mask[in.second]=true;
    }
    vector<float> trunkfreq;
    for(int i=0;i<mask.size();i++)
        if(mask[i])
            trunkfreq.push_back(m_freqs[i]);
    sort(trunkfreq.begin(),trunkfreq.end());
    for(int i=0;i<trunkfreq.size();i++)
        cout<<"Inlier: "<<trunkfreq[i]<<endl;
}

vector<float> FrequencySieve::getOutlierFrequencies()
{
    vector<bool> mask(m_freqs.size(),false);
    vector<int > inl=inliers();
    for(int i=0;i<inl.size();i++)
    {
        pair<int,int> in=m_idxs[inl[i]];
        mask[in.first]=true;
        mask[in.second]=true;
    }
    vector<float> trunkfreq;
    for(int i=0;i<mask.size();i++)
        if(!mask[i])
            trunkfreq.push_back(m_freqs[i]);
    return trunkfreq;
}

vector<float> FrequencySieve::getOutlierAmplitudes()
{
    vector<bool> mask(m_freqs.size(),false);
    vector<int > inl=inliers();
    for(int i=0;i<inl.size();i++)
    {
        pair<int,int> in=m_idxs[inl[i]];
        mask[in.first]=true;
        mask[in.second]=true;
    }
    vector<float> trunkamp;
    for(int i=0;i<mask.size();i++)
        if(!mask[i])
            trunkamp.push_back(m_amplitudes[i]);
    return trunkamp;
}
