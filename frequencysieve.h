#ifndef FREQUENCYSIEVE_H
#define FREQUENCYSIEVE_H
#include <multiRansac.hxx>
#include <vector>
#include <algorithm>
using namespace std;
class FrequencySieve: public Ransac<pair<int,int>, float>
{
public:
    FrequencySieve(vector<float >freqs,vector<float > amplitudes,double thresh);


    virtual vector<double> errorfunction(pair<int,int> dataid,const vector<float> &currF0,const vector<int>& inliers);

    virtual vector<float> computeParameters(const vector<pair<int,int> > &data,int numunknowns);
    float computedF0()
    {
        return unknowns()[0];
    }
    void printInliers();

    vector<float> getOutlierFrequencies();

    vector<float> getOutlierAmplitudes();

private:
    vector<float> m_freqs;
    vector<pair<int,int> > m_idxs;
    vector<float> m_amplitudes;
};

#endif // FREQUENCYSIEVE_H
