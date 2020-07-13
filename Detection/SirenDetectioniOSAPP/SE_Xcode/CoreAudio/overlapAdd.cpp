////
///
///
//
//  overlapAdd.cpp
//  frameOverlap
//
//  Created by Shankar, Nikhil on 8/28/19.
//  Copyright © 2019 default. All rights reserved.
//

#include "overlapAdd.hpp"
#include <math.h>
#include <algorithm>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define M_D 0.6100
#define M_V 0.7050
#define FRAMESEARCH 10
#define WINDOWSEARCH 20
#define MAX_VALUE 1E+09
#define TRAININGFRAMES 6
#define vadThreshold 0.15
#define BETA 0.85
using namespace std;
void overlapAdd::init()
{
    in_prev = new float [frameS];
    outputOld = new float [frameS];
    win = new float [frameS*2];
    de = new float [(frameS*2)-1];
    dewin = new float [frameS*2];
    invWin = new float [frameS];
    inputBuffer = new float [frameS*2];
    xReal1 = new float [nFFT];
    xImag1 = new float [nFFT];
    yReal = new float [nFFT];
    yImag = new float [nFFT];
    xReal2 = new float [nFFT];
    xImag2 = new float [nFFT];
    cosine = new float [nFFT/2];
    sine = new float [nFFT];
    mag = new float [nFFT/2 +1];
    log_mag = new float [nFFT/2 + 1];
    noiseEstimatePrev = new float [nFFT/2 +1];
    Pprev = new float [nFFT/2 +1];
    meanP = new float [nFFT/2 +1];
    meanP2 = new float [nFFT/2 +1];
    alphaCPrev = 0;
    alphaHat = new float [nFFT/2 +1];
    P = new float [nFFT/2 +1];
    beta = new float [nFFT/2 +1];
    varP = new float [nFFT/2 +1];
    invQeq = new float [nFFT/2 +1];
    tilda_QD = new float [nFFT/2 +1];
    tilda_QV = new float [nFFT/2 +1];
    Bmin = new float [nFFT/2 +1];
    Bminsub = new float [nFFT/2 +1];
    PBC = new float [nFFT/2 +1];
    kMode = new int [nFFT/2 +1];
    actmin = new float [nFFT/2 +1];
    actminSub = new float [nFFT/2 +1];
    limflag = new int [nFFT/2 +1];
    Pmin = new float [nFFT/2 +1];
    p = new int [nFFT/2 +1];
    temp1 = new int [nFFT/2 +1];
    temp2 = new int [nFFT/2 +1];
    noiseEstimate = new float [nFFT/2 +1];
    noiseMean = new float [nFFT/2 +1];
    noisePower = new float [nFFT/2 +1];
    noiseMu = new float [nFFT/2 +1];
    postSNR = new float [nFFT/2 +1];
    aprioriSNR = new float [nFFT/2 +1];
    prevOutput = new float [nFFT/2 +1];
    postSNR1 = new float [nFFT/2 +1];
    aprioriSNR1 = new float [nFFT/2 +1];
    gain = new float [nFFT/2 +1];
    prevGain = new float [nFFT/2 +1];
    ensig = new float [nFFT/2 +1];
    logSigma = new float [nFFT/2 +1];
    finalSig = new float [nFFT];
    angle = new float [nFFT];
    
    ///siren detect parameters
    log_mag = new float [nFFT/2 + 1];

    for(int i=0;i<frameS;++i)
    {
        in_prev[i]=0;
        outputOld[i]=0;
    }
    for(int i=0;i<nFFT/2;++i)
    {
        cosine[i]=0;
        sine[i]=0;
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        actmin[i] = MAX_VALUE;
        actminSub[i] = MAX_VALUE;
        kMode[i] = 0;
        limflag[i] = 0;
        noiseMean[i] = 0;
        prevGain[i] = 0;
    }
    subwc = WINDOWSEARCH;
    count1 =0;
    count0 =0;
}

overlapAdd::~overlapAdd()
{
    delete [] in_prev;
    delete [] inputBuffer;
    delete [] win;
    delete [] xReal1;
    delete [] xImag1;
    delete [] xReal2;
    delete [] xImag2;
    delete [] yReal;
    delete [] yImag;
    delete [] outputOld;
    delete [] invWin;
    delete [] mag;
    delete [] noiseEstimatePrev;
    delete [] Pprev;
    delete [] meanP;
    delete [] meanP2;
    delete [] alphaHat;
    delete [] P;
    delete [] beta;
    delete [] varP;
    delete [] invQeq;
    delete [] tilda_QV;
    delete [] tilda_QD;
    delete [] Bmin;
    delete [] Bminsub;
    delete [] PBC;
    delete [] kMode;
    delete [] actmin;
    delete [] actminSub;
    delete [] limflag;
    delete [] Pmin;
    delete [] p;
    delete [] temp1;
    delete [] temp2;
    delete [] noiseEstimate;
    delete [] noiseMean;
    delete [] noisePower;
    delete [] noiseMu;
    delete [] postSNR;
    delete [] aprioriSNR;
    delete [] prevOutput;
    delete [] postSNR1;
    delete [] aprioriSNR1;
    delete [] gain;
    delete [] prevGain;
    delete [] ensig;
    delete [] logSigma;
    delete [] finalSig;
    delete [] angle;
    delete [] de;
    delete [] dewin;
}

void overlapAdd::process(const float *input, float *output, const int frameCount)
{
    float sumYenergy = 0;
    float sumPenergy = 0;
    float sumNoise = 0;
    float sumInv = 0;
    float sumLogSigma = 0;
    if(frameCount==0)
    {
        for (int i = 0; i < (frameS*2); i++)
        {
            win[i] = 0.54 - 0.46 * (cos(2 * M_PI*(i) / ((frameS*2)-1)));
        }
        for (int i = 0; i < (frameS*2)-1; i++)
        {
            de[i] = 0.5 * (1 - cosf(2 * M_PI*(i + 1) / (((frameS*2)-1) + 1)));
        }
        for (int i = 0; i < frameS; i++)
        {
            dewin[i] = de[i]/win[i];
        }
        for (int i = frameS; i < (frameS*2); i++)
        {
            dewin[i] = de[i-1]/win[i];
        }
    }
    for(int i=0;i<frameS;++i)
    {
        inputBuffer[i] = in_prev[i];
        in_prev[i]=input[i];
        inputBuffer[i+frameS] = input[i];
    }
    for(int i=0;i<frameS*2;++i)
    {
        inputBuffer[i] = inputBuffer[i]*win[i];
    }
    FFT(inputBuffer,xReal1,xImag1,nFFT);
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        mag[i] = (xReal1[i] * xReal1[i]) + (xImag1[i] * xImag1[i]);
    }
    if(frameCount==0)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            noiseEstimatePrev[i] = mag[i];
            Pprev[i] = mag[i];
            meanP[i] = mag[i];
            meanP2[i] = pow(mag[i],2);
        }
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        sumYenergy += mag[i];
        sumPenergy += Pprev[i];
    }
    if(sumYenergy==0)
    {
        sumYenergy = 1E-09;
    }
    alphaC = 1/(1+pow((sumPenergy/sumYenergy - 1),2));
    alphaC = (0.7 * alphaCPrev)+(0.3 * MAX(alphaC,0.3));
    alphaCPrev = alphaC;
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        alphaHat[i] = (0.96 * alphaC/(1 + pow((Pprev[i]/noiseEstimatePrev[i]-1),2)));
        sumNoise += noiseEstimatePrev[i];
    }
    snr = sumPenergy/sumNoise;
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        alphaHat[i] = MAX(alphaHat[i], MIN(0.3, pow(snr,(nFFT/(2*0.064*fs)))));
        P[i] = (alphaHat[i] * Pprev[i]) + ((1 - alphaHat[i])*mag[i]);
        beta[i] = MIN(pow(alphaHat[i],2), 0.8);
        meanP[i] = (beta[i] * meanP[i]) + ((1-beta[i])*P[i]);
        meanP2[i] = (beta[i] * meanP2[i]) + ((1-beta[i])*pow(P[i],2));
        varP[i] = meanP2[i] - pow(meanP[i],2);
        invQeq[i] = varP[i]/(2 * pow(noiseEstimatePrev[i],2));
        invQeq[i] = MAX(MIN(invQeq[i], 0.5), 1/(20.0*(frameCount+1)));
        tilda_QD[i] = (1/invQeq[i] - 2*M_D)/(1-M_D);
        tilda_QV[i] = (1/invQeq[i] - 2*M_V)/(1-M_V);
        Bmin[i] = 1 + 2*(FRAMESEARCH-1)/tilda_QD[i];
        Bminsub[i] = 1 + 2*(WINDOWSEARCH-1)/tilda_QV[i];
        sumInv += invQeq[i];
    }
    meanInv = sumInv/(nFFT/2 +1);
    Bc = 1 + 2.12*sqrt(meanInv);
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        PBC[i] = P[i] * Bmin[i] * Bc;
        if(PBC[i] < actmin[i])
        {
            kMode[i] = 1;
            actmin[i] = PBC[i];
            actminSub[i] = P[i] * Bminsub[i] * Bc;
        }
        else
        {
            kMode[i] = 0;
        }
    }
    if(subwc == WINDOWSEARCH)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            if(kMode[i]==1)
            {
                limflag[i] = 0;
            }
            Pmin[i] = actmin[i];
            
        }
        if(meanInv<0.03)
        {
            noiseSlopeMax = 8;
        }
        else if ((meanInv>=0.03) && (meanInv<0.05))
        {
            noiseSlopeMax = 4;
        }
        else if ((meanInv>=0.05) && (meanInv<0.06))
        {
            noiseSlopeMax = 2;
        }
        else if (meanInv>=0.06)
        {
            noiseSlopeMax = 1.2;
        }
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            if(actminSub[i] < noiseSlopeMax*Pmin[i])
            {
                temp1[i] = 1;
            }
            else
            {
                temp1[i] = 0;
            }
            if(actminSub[i] > Pmin[i])
            {
                temp2[i] = 1;
            }
            else
            {
                temp2[i] = 0;
            }
            p[i] = (limflag[i] & temp1[i]) & temp2[i];
            if(p[i] == 1)
            {
                Pmin[i] = actminSub[i];
                actmin[i] = actminSub[i];
            }
        }
        subwc = 1;
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            actmin[i] = MAX_VALUE;
            actminSub[i] = MAX_VALUE;
            noiseEstimate[i] = Pmin[i];
        }
    }
    else
    {
        if(subwc > 1)
        {
            for(int i=0;i<nFFT/2 + 1;++i)
            {
                if(kMode[i]==1)
                {
                    limflag[i] = 1;
                }
                noiseEstimate[i] = MIN(actminSub[i],Pmin[i]);
                Pmin[i] = noiseEstimate[i];
            }
        }
        else
        {
            for(int i=0;i<nFFT/2 + 1;++i)
            {
                noiseEstimate[i] = MIN(actminSub[i],Pmin[i]);
            }
        }
        subwc = subwc + 1;
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        Pprev[i] = P[i];
        noiseEstimatePrev[i] = noiseEstimate[i];
    }
    
    if(frameCount<TRAININGFRAMES)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            noiseMean[i] += sqrt(mag[i]);
            noisePower[i] = mag[i];
        }
    }
    if(frameCount == TRAININGFRAMES-1)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            noiseMu[i] = noiseMean[i]/TRAININGFRAMES;
            noisePower[i] = pow(noiseMu[i],2);
        }
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        postSNR[i] = P[i]/noiseEstimate[i];
    }
    if(frameCount==0)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            aprioriSNR[i] = (0.98) + (1-0.98)*MAX(postSNR[i]-1, 0);
        }
    }
    else
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            aprioriSNR[i] = (0.98 * (prevOutput[i]/noiseEstimate[i]))+ (1-0.98)*MAX(postSNR[i]-1, 0);
            aprioriSNR[i] = MAX((float)pow(10,((float)-25 /(float)10)), aprioriSNR[i]);
        }
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        logSigma[i] = postSNR[i] * aprioriSNR[i] / (1 + aprioriSNR[i]) - log(1 + aprioriSNR[i]);
        sumLogSigma += logSigma[i];
    }
    vadDecision = sumLogSigma/(nFFT/2 +1);
    prevVad = vad;
    if(vadDecision<vadThreshold)
    {
        vad = 0;
    }
    else
    {
        vad = 1;
    }
    if(frameCount>0)
    {
        if(vad == 1)
        {
            count1 = 0;
            ++count0;
            if(count0<1)
            {
                vad = prevVad;
            }
            else
            {
                vad = 1;
            }
        }
        else if (vad == 0)
        {
            count0=0;
            ++count1;
            if(count1<15)
            {
                vad = prevVad;
            }
            else
            {
                vad = 0;
            }
        }
    }
    if((vad==0) && (frameCount>TRAININGFRAMES-1))
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            noisePower[i] = (0.99 * noisePower[i]) + ((1-0.99) * noiseEstimate[i]);
        }
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        postSNR1[i] = MIN((mag[i]/noisePower[i]),40);
    }
    if(frameCount==0)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            aprioriSNR1[i] = (0.98) + (1-0.98)*MAX(postSNR1[i]-1, 0);
        }
    }
    else
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            aprioriSNR1[i] = (0.98 * (prevOutput[i]/noisePower[i]))+ (1-0.98)*MAX(postSNR1[i]-1, 0);
            aprioriSNR1[i] = MAX((float)pow(10,((float)-25 /(float)10)), aprioriSNR1[i]);
        }
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        gain[i] = sqrt(aprioriSNR1[i])/(BETA + sqrt(aprioriSNR1[i]));
    }
    float minGain = gain[0];
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        if(minGain>gain[i])
        {
            minGain=gain[i];
        }
    }
    if(vad==0)
    {
        for(int i=0;i<nFFT/2 + 1;++i)
        {
            if(gain[i] > (minGain * 4))
            {
                gain[i] = gain[i]/1.5;
            }
        }
    }
    for(int i=44;i<nFFT/2 + 1;++i)
    {
        gain[i] = gain[i] / 1.05;
    }
    for(int i=0;i<20;++i)
    {
        gain[i] = gain[i] * 1.1;
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        gain[i] = ((1-0.5)*gain[i]) + (0.5 * prevGain[i]);
        prevGain[i] = gain[i];
        ensig[i] = sqrt(mag[i]) * gain[i];
        prevOutput[i] = pow(ensig[i],2);
    }
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        finalSig[i] = ensig[i];
    }
    std::reverse(&ensig[0],&ensig[nFFT/2 + 1]);
    int j=1;
    for(int i=nFFT/2 + 1;i<nFFT;++i)
    {
        finalSig[i] = ensig[j];
        ++j;
    }
    for(int i=0;i<nFFT;++i)
    {
        angle[i] = atan2(xImag1[i], xReal1[i]);
        xReal1[i] = cos(angle[i])*finalSig[i];
        xImag1[i] = sin(angle[i])*finalSig[i];
    }
    IFFT(xReal1, xImag1, yReal, yImag, nFFT);
    for(int i=0;i<nFFT;++i)
    {
        yReal[i] = yReal[i] * dewin[i];
    }
    for(int i=0; i<frameS; i++)
    {
        output[i]= outputOld[i]+yReal[i];
        outputOld[i]=yReal[i+nFFT/2];
    }
}

void overlapAdd::siren_det(const float *input, float *siren_fea, const int frameCount)
{

    if(frameCount==0)
    {
        for (int i = 0; i < (frameS*2); i++)
        {
            win[i] = 0.54 - 0.46 * (cos(2 * M_PI*(i) / ((frameS*2)-1)));
        }
        for (int i = 0; i < (frameS*2)-1; i++)
        {
            de[i] = 0.5 * (1 - cosf(2 * M_PI*(i + 1) / (((frameS*2)-1) + 1)));
        }
        for (int i = 0; i < frameS; i++)
        {
            dewin[i] = de[i]/win[i];
        }
        for (int i = frameS; i < (frameS*2); i++)
        {
            dewin[i] = de[i-1]/win[i];
        }
    }
    for(int i=0;i<frameS;++i)
    {
        inputBuffer[i] = in_prev[i];
        in_prev[i]=input[i];
        inputBuffer[i+frameS] = input[i];
    }
    for(int i=0;i<frameS*2;++i)
    {
        inputBuffer[i] = inputBuffer[i]*win[i];
    }
    FFT(inputBuffer,xReal1,xImag1,nFFT);
    for(int i=0;i<nFFT/2 + 1;++i)
    {
        log_mag[i] = 20 * log10((xReal1[i] * xReal1[i]) + (xImag1[i] * xImag1[i]));
    }

    for(int i=0; i<frameS; i++)
    {
        siren_fea[i]= log_mag[i];
    }
}



void overlapAdd::FFT(const float* input, float *outputReal, float *outputImag, const int nFFT)
{
    float arg;
    for (int i = 0; i<nFFT / 2; i++)
    {
        arg = -2 * M_PI*i / nFFT;
        cosine[i] = cos(arg);
        sine[i] = sin(arg);
    }
    int i, j, k, L, m, n, o, p, q;
    float tempReal, tempImaginary, cos, sin, xt, yt;
    k = nFFT;
    for (i = 0; i<k; i++)
    {
        outputReal[i] = input[i];
        outputImag[i] = 0;
    }
    
    j = 0;
    m = k / 2;
    //bit reversal
    for (i = 1; i<(k - 1); i++)
    {
        L = m;
        //L = pow(2,ceil(log2(m)));
        while (j >= L)
        {
            j = j - L;
            L = L / 2;
        }
        j = j + L;
        if (i<j)
        {
            tempReal = outputReal[i];
            tempImaginary = outputImag[i];
            outputReal[i] = outputReal[j];
            outputImag[i] = outputImag[j];
            outputReal[j] = tempReal;
            outputImag[j] = tempImaginary;
        }
    }
    L = 0;
    m = 1;
    n = k / 2;
    //computation
    for (i = k; i>1; i = (i >> 1))
    {
        L = m;
        m = 2 * m;
        o = 0;
        for (j = 0; j<L; j++)
        {
            cos = cosine[o];
            sin = sine[o];
            o = o + n;
            for (p = j; p<k; p = p + m)
            {
                q = p + L;
                xt = cos*outputReal[q] - sin*outputImag[q];
                yt = sin*outputReal[q] + cos*outputImag[q];
                outputReal[q] = (outputReal[p] - xt);
                outputImag[q] = (outputImag[p] - yt);
                outputReal[p] = (outputReal[p] + xt);
                outputImag[p] = (outputImag[p] + yt);
            }
        }
        n = n >> 1;
    }
}

void overlapAdd::IFFT(const float* inputReal, const float* inputImag, float *outputReal, float *outputImag, const int nFFT)
{
    float arg;
    for (int i = 0; i<nFFT / 2; i++)
    {
        arg = -2 * M_PI*i / nFFT;
        cosine[i] = cos(arg);
        sine[i] = sin(arg);
    }
    int i, j, k, L, m, n, o, p, q;
    float tempReal, tempImaginary, cos, sin, xt, yt;
    k = nFFT;
    for (i = 0; i<k; i++)
    {
        outputReal[i] = inputReal[i];
        outputImag[i] = (-1)*inputImag[i];
    }
    
    j = 0;
    m = k / 2;
    //bit reversal
    for (i = 1; i<(k - 1); i++)
    {
        L = m;
        while (j >= L)
        {
            j = j - L;
            L = L / 2;
        }
        j = j + L;
        if (i<j)
        {
            tempReal = outputReal[i];
            tempImaginary = outputImag[i];
            outputReal[i] = outputReal[j];
            outputImag[i] = outputImag[j];
            outputReal[j] = tempReal;
            outputImag[j] = tempImaginary;
        }
    }
    L = 0;
    m = 1;
    n = k / 2;
    //computation
    for (i = k; i>1; i = (i >> 1))
    {
        L = m;
        m = 2 * m;
        o = 0;
        for (j = 0; j<L; j++)
        {
            cos = cosine[o];
            sin = sine[o];
            o = o + n;
            for (p = j; p<k; p = p + m)
            {
                q = p + L;
                xt = cos*outputReal[q] - sin*outputImag[q];
                yt = sin*outputReal[q] + cos*outputImag[q];
                outputReal[q] = (outputReal[p] - xt);
                outputImag[q] = (outputImag[p] - yt);
                outputReal[p] = (outputReal[p] + xt);
                outputImag[p] = (outputImag[p] + yt);
            }
        }
        n = n >> 1;
    }
    for (i = 0; i<k; i++)
    {
        outputReal[i] = outputReal[i] / k;
        outputImag[i] = outputImag[i] / k;
    }
}

////  overlapAdd.cpp
////  frameOverlap
////
////  Created by Shankar, Nikhil on 8/28/19.
////  Copyright © 2019 default. All rights reserved.
////
//
//#include "overlapAdd.hpp"
//#include <math.h>
//#include <algorithm>
//
//#define MIN(a,b) (((a)<(b))?(a):(b))
//#define MAX(a,b) (((a)>(b))?(a):(b))
//#define M_D 0.6100
//#define M_V 0.7050
//#define FRAMESEARCH 10
//#define WINDOWSEARCH 20
//#define MAX_VALUE 1E+09
//#define TRAININGFRAMES 6
//#define vadThreshold 0.15
//#define BETA 0.85
//using namespace std;
//void overlapAdd::init()
//{
//    in_prev = new float [frameS];
//    outputOld = new float [frameS];
//    win = new float [frameS*2];
//    de = new float [(frameS*2)-1];
//    dewin = new float [frameS*2];
//    invWin = new float [frameS];
//    inputBuffer = new float [frameS*2];
//    xReal1 = new float [nFFT];
//    xImag1 = new float [nFFT];
//    yReal = new float [nFFT];
//    yImag = new float [nFFT];
//    xReal2 = new float [nFFT];
//    xImag2 = new float [nFFT];
//    cosine = new float [nFFT/2];
//    sine = new float [nFFT];
//    mag = new float [nFFT/2 +1];
//    noiseEstimatePrev = new float [nFFT/2 +1];
//    Pprev = new float [nFFT/2 +1];
//    meanP = new float [nFFT/2 +1];
//    meanP2 = new float [nFFT/2 +1];
//    alphaCPrev = 0;
//    alphaHat = new float [nFFT/2 +1];
//    P = new float [nFFT/2 +1];
//    beta = new float [nFFT/2 +1];
//    varP = new float [nFFT/2 +1];
//    invQeq = new float [nFFT/2 +1];
//    tilda_QD = new float [nFFT/2 +1];
//    tilda_QV = new float [nFFT/2 +1];
//    Bmin = new float [nFFT/2 +1];
//    Bminsub = new float [nFFT/2 +1];
//    PBC = new float [nFFT/2 +1];
//    kMode = new int [nFFT/2 +1];
//    actmin = new float [nFFT/2 +1];
//    actminSub = new float [nFFT/2 +1];
//    limflag = new int [nFFT/2 +1];
//    Pmin = new float [nFFT/2 +1];
//    p = new int [nFFT/2 +1];
//    temp1 = new int [nFFT/2 +1];
//    temp2 = new int [nFFT/2 +1];
//    noiseEstimate = new float [nFFT/2 +1];
//    noiseMean = new float [nFFT/2 +1];
//    noisePower = new float [nFFT/2 +1];
//    noiseMu = new float [nFFT/2 +1];
//    postSNR = new float [nFFT/2 +1];
//    aprioriSNR = new float [nFFT/2 +1];
//    prevOutput = new float [nFFT/2 +1];
//    postSNR1 = new float [nFFT/2 +1];
//    aprioriSNR1 = new float [nFFT/2 +1];
//    gain = new float [nFFT/2 +1];
//    prevGain = new float [nFFT/2 +1];
//    ensig = new float [nFFT/2 +1];
//    logSigma = new float [nFFT/2 +1];
//    finalSig = new float [nFFT];
//    angle = new float [nFFT];
//
//    for(int i=0;i<frameS;++i)
//    {
//        in_prev[i]=0;
//        outputOld[i]=0;
//    }
//    for(int i=0;i<nFFT/2;++i)
//    {
//        cosine[i]=0;
//        sine[i]=0;
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        actmin[i] = MAX_VALUE;
//        actminSub[i] = MAX_VALUE;
//        kMode[i] = 0;
//        limflag[i] = 0;
//        noiseMean[i] = 0;
//        prevGain[i] = 0;
//    }
//    subwc = WINDOWSEARCH;
//    count1 =0;
//    count0 =0;
//}
//
//overlapAdd::~overlapAdd()
//{
//    delete [] in_prev;
//    delete [] inputBuffer;
//    delete [] win;
//    delete [] xReal1;
//    delete [] xImag1;
//    delete [] xReal2;
//    delete [] xImag2;
//    delete [] yReal;
//    delete [] yImag;
//    delete [] outputOld;
//    delete [] invWin;
//    delete [] mag;
//    delete [] noiseEstimatePrev;
//    delete [] Pprev;
//    delete [] meanP;
//    delete [] meanP2;
//    delete [] alphaHat;
//    delete [] P;
//    delete [] beta;
//    delete [] varP;
//    delete [] invQeq;
//    delete [] tilda_QV;
//    delete [] tilda_QD;
//    delete [] Bmin;
//    delete [] Bminsub;
//    delete [] PBC;
//    delete [] kMode;
//    delete [] actmin;
//    delete [] actminSub;
//    delete [] limflag;
//    delete [] Pmin;
//    delete [] p;
//    delete [] temp1;
//    delete [] temp2;
//    delete [] noiseEstimate;
//    delete [] noiseMean;
//    delete [] noisePower;
//    delete [] noiseMu;
//    delete [] postSNR;
//    delete [] aprioriSNR;
//    delete [] prevOutput;
//    delete [] postSNR1;
//    delete [] aprioriSNR1;
//    delete [] gain;
//    delete [] prevGain;
//    delete [] ensig;
//    delete [] logSigma;
//    delete [] finalSig;
//    delete [] angle;
//    delete [] de;
//    delete [] dewin;
//}
//
//void overlapAdd::process(const float *input, float *output, const int frameCount)
//{
//    float sumYenergy = 0;
//    float sumPenergy = 0;
//    float sumNoise = 0;
//    float sumInv = 0;
//    float sumLogSigma = 0;
//    if(frameCount==0)
//    {
//        for (int i = 0; i < (frameS*2); i++)
//        {
//            win[i] = 0.54 - 0.46 * (cos(2 * M_PI*(i) / ((frameS*2)-1)));
//        }
//        for (int i = 0; i < (frameS*2)-1; i++)
//        {
//            de[i] = 0.5 * (1 - cosf(2 * M_PI*(i + 1) / (((frameS*2)-1) + 1)));
//        }
//        for (int i = 0; i < frameS; i++)
//        {
//            dewin[i] = de[i]/win[i];
//        }
//        for (int i = frameS; i < (frameS*2); i++)
//        {
//            dewin[i] = de[i-1]/win[i];
//        }
//    }
//    for(int i=0;i<frameS;++i)
//    {
//        inputBuffer[i] = in_prev[i];
//        in_prev[i]=input[i];
//        inputBuffer[i+frameS] = input[i];
//    }
//    for(int i=0;i<frameS*2;++i)
//    {
//        inputBuffer[i] = inputBuffer[i]*win[i];
//    }
//    FFT(inputBuffer,xReal1,xImag1,nFFT);
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
////        mag[i] = 10*log10 ((xReal1[i] * xReal1[i]) + (xImag1[i] * xImag1[i]));
//        mag[i] = ((xReal1[i] * xReal1[i]) + (xImag1[i] * xImag1[i]));
//
//    }
//    //SE
//    if(frameCount==0)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            noiseEstimatePrev[i] = mag[i];
//            Pprev[i] = mag[i];
//            meanP[i] = mag[i];
//            meanP2[i] = pow(mag[i],2);
//        }
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        sumYenergy += mag[i];
//        sumPenergy += Pprev[i];
//    }
//    if(sumYenergy==0)
//    {
//        sumYenergy = 1E-09;
//    }
//    alphaC = 1/(1+pow((sumPenergy/sumYenergy - 1),2));
//    alphaC = (0.7 * alphaCPrev)+(0.3 * MAX(alphaC,0.3));
//    alphaCPrev = alphaC;
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        alphaHat[i] = (0.96 * alphaC/(1 + pow((Pprev[i]/noiseEstimatePrev[i]-1),2)));
//        sumNoise += noiseEstimatePrev[i];
//    }
//    snr = sumPenergy/sumNoise;
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        alphaHat[i] = MAX(alphaHat[i], MIN(0.3, pow(snr,(nFFT/(2*0.064*fs)))));
//        P[i] = (alphaHat[i] * Pprev[i]) + ((1 - alphaHat[i])*mag[i]);
//        beta[i] = MIN(pow(alphaHat[i],2), 0.8);
//        meanP[i] = (beta[i] * meanP[i]) + ((1-beta[i])*P[i]);
//        meanP2[i] = (beta[i] * meanP2[i]) + ((1-beta[i])*pow(P[i],2));
//        varP[i] = meanP2[i] - pow(meanP[i],2);
//        invQeq[i] = varP[i]/(2 * pow(noiseEstimatePrev[i],2));
//        invQeq[i] = MAX(MIN(invQeq[i], 0.5), 1/(20.0*(frameCount+1)));
//        tilda_QD[i] = (1/invQeq[i] - 2*M_D)/(1-M_D);
//        tilda_QV[i] = (1/invQeq[i] - 2*M_V)/(1-M_V);
//        Bmin[i] = 1 + 2*(FRAMESEARCH-1)/tilda_QD[i];
//        Bminsub[i] = 1 + 2*(WINDOWSEARCH-1)/tilda_QV[i];
//        sumInv += invQeq[i];
//    }
//    meanInv = sumInv/(nFFT/2 +1);
//    Bc = 1 + 2.12*sqrt(meanInv);
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        PBC[i] = P[i] * Bmin[i] * Bc;
//        if(PBC[i] < actmin[i])
//        {
//            kMode[i] = 1;
//            actmin[i] = PBC[i];
//            actminSub[i] = P[i] * Bminsub[i] * Bc;
//        }
//        else
//        {
//            kMode[i] = 0;
//        }
//    }
//    if(subwc == WINDOWSEARCH)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            if(kMode[i]==1)
//            {
//                limflag[i] = 0;
//            }
//            Pmin[i] = actmin[i];
//
//        }
//        if(meanInv<0.03)
//        {
//            noiseSlopeMax = 8;
//        }
//        else if ((meanInv>=0.03) && (meanInv<0.05))
//        {
//            noiseSlopeMax = 4;
//        }
//        else if ((meanInv>=0.05) && (meanInv<0.06))
//        {
//            noiseSlopeMax = 2;
//        }
//        else if (meanInv>=0.06)
//        {
//            noiseSlopeMax = 1.2;
//        }
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            if(actminSub[i] < noiseSlopeMax*Pmin[i])
//            {
//                temp1[i] = 1;
//            }
//            else
//            {
//                temp1[i] = 0;
//            }
//            if(actminSub[i] > Pmin[i])
//            {
//                temp2[i] = 1;
//            }
//            else
//            {
//                temp2[i] = 0;
//            }
//            p[i] = (limflag[i] & temp1[i]) & temp2[i];
//            if(p[i] == 1)
//            {
//                Pmin[i] = actminSub[i];
//                actmin[i] = actminSub[i];
//            }
//        }
//        subwc = 1;
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            actmin[i] = MAX_VALUE;
//            actminSub[i] = MAX_VALUE;
//            noiseEstimate[i] = Pmin[i];
//        }
//    }
//    else
//    {
//        if(subwc > 1)
//        {
//            for(int i=0;i<nFFT/2 + 1;++i)
//            {
//                if(kMode[i]==1)
//                {
//                    limflag[i] = 1;
//                }
//                noiseEstimate[i] = MIN(actminSub[i],Pmin[i]);
//                Pmin[i] = noiseEstimate[i];
//            }
//        }
//        else
//        {
//            for(int i=0;i<nFFT/2 + 1;++i)
//            {
//                noiseEstimate[i] = MIN(actminSub[i],Pmin[i]);
//            }
//        }
//        subwc = subwc + 1;
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        Pprev[i] = P[i];
//        noiseEstimatePrev[i] = noiseEstimate[i];
//    }
//
//    if(frameCount<TRAININGFRAMES)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            noiseMean[i] += sqrt(mag[i]);
//            noisePower[i] = mag[i];
//        }
//    }
//    if(frameCount == TRAININGFRAMES-1)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            noiseMu[i] = noiseMean[i]/TRAININGFRAMES;
//            noisePower[i] = pow(noiseMu[i],2);
//        }
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        postSNR[i] = P[i]/noiseEstimate[i];
//    }
//    if(frameCount==0)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            aprioriSNR[i] = (0.98) + (1-0.98)*MAX(postSNR[i]-1, 0);
//        }
//    }
//    else
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            aprioriSNR[i] = (0.98 * (prevOutput[i]/noiseEstimate[i]))+ (1-0.98)*MAX(postSNR[i]-1, 0);
//            aprioriSNR[i] = MAX((float)pow(10,((float)-25 /(float)10)), aprioriSNR[i]);
//        }
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        logSigma[i] = postSNR[i] * aprioriSNR[i] / (1 + aprioriSNR[i]) - log(1 + aprioriSNR[i]);
//        sumLogSigma += logSigma[i];
//    }
//    vadDecision = sumLogSigma/(nFFT/2 +1);
//    prevVad = vad;
//    if(vadDecision<vadThreshold)
//    {
//        vad = 0;
//    }
//    else
//    {
//        vad = 1;
//    }
//    if(frameCount>0)
//    {
//        if(vad == 1)
//        {
//            count1 = 0;
//            ++count0;
//            if(count0<1)
//            {
//                vad = prevVad;
//            }
//            else
//            {
//                vad = 1;
//            }
//        }
//        else if (vad == 0)
//        {
//            count0=0;
//            ++count1;
//            if(count1<15)
//            {
//                vad = prevVad;
//            }
//            else
//            {
//                vad = 0;
//            }
//        }
//    }
//    if((vad==0) && (frameCount>TRAININGFRAMES-1))
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            noisePower[i] = (0.99 * noisePower[i]) + ((1-0.99) * noiseEstimate[i]);
//        }
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        postSNR1[i] = MIN((mag[i]/noisePower[i]),40);
//    }
//    if(frameCount==0)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            aprioriSNR1[i] = (0.98) + (1-0.98)*MAX(postSNR1[i]-1, 0);
//        }
//    }
//    else
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            aprioriSNR1[i] = (0.98 * (prevOutput[i]/noisePower[i]))+ (1-0.98)*MAX(postSNR1[i]-1, 0);
//            aprioriSNR1[i] = MAX((float)pow(10,((float)-25 /(float)10)), aprioriSNR1[i]);
//        }
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        gain[i] = sqrt(aprioriSNR1[i])/(BETA + sqrt(aprioriSNR1[i]));
//    }
//    float minGain = gain[0];
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        if(minGain>gain[i])
//        {
//            minGain=gain[i];
//        }
//    }
//    if(vad==0)
//    {
//        for(int i=0;i<nFFT/2 + 1;++i)
//        {
//            if(gain[i] > (minGain * 4))
//            {
//                gain[i] = gain[i]/1.5;
//            }
//        }
//    }
//    for(int i=44;i<nFFT/2 + 1;++i)
//    {
//        gain[i] = gain[i] / 1.05;
//    }
//    for(int i=0;i<20;++i)
//    {
//        gain[i] = gain[i] * 1.1;
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        gain[i] = ((1-0.5)*gain[i]) + (0.5 * prevGain[i]);
//        prevGain[i] = gain[i];
//        ensig[i] = sqrt(mag[i]) * gain[i];
//        prevOutput[i] = pow(ensig[i],2);
//    }
//    for(int i=0;i<nFFT/2 + 1;++i)
//    {
//        finalSig[i] = ensig[i];
//    }
//    std::reverse(&ensig[0],&ensig[nFFT/2 + 1]);
//    int j=1;
//    for(int i=nFFT/2 + 1;i<nFFT;++i)
//    {
//        finalSig[i] = ensig[j];
//        ++j;
//    }
//    for(int i=0;i<nFFT;++i)
//    {
//        angle[i] = atan2(xImag1[i], xReal1[i]);
//        xReal1[i] = cos(angle[i])*finalSig[i];
//        xImag1[i] = sin(angle[i])*finalSig[i];
//    }
//    IFFT(xReal1, xImag1, yReal, yImag, nFFT);
//    for(int i=0;i<nFFT;++i)
//    {
//        yReal[i] = yReal[i] * dewin[i];
//    }
//    for(int i=0; i<frameS; i++)
//    {
//        output[i]= outputOld[i]+yReal[i];
//        outputOld[i]=yReal[i+nFFT/2];
////        output[i]=mag[i];
//    }
//}
//
//void overlapAdd::FFT(const float* input, float *outputReal, float *outputImag, const int nFFT)
//{
//    float arg;
//    for (int i = 0; i<nFFT / 2; i++)
//    {
//        arg = -2 * M_PI*i / nFFT;
//        cosine[i] = cos(arg);
//        sine[i] = sin(arg);
//    }
//    int i, j, k, L, m, n, o, p, q;
//    float tempReal, tempImaginary, cos, sin, xt, yt;
//    k = nFFT;
//    for (i = 0; i<k; i++)
//    {
//        outputReal[i] = input[i];
//        outputImag[i] = 0;
//    }
//
//    j = 0;
//    m = k / 2;
//    //bit reversal
//    for (i = 1; i<(k - 1); i++)
//    {
//        L = m;
//        //L = pow(2,ceil(log2(m)));
//        while (j >= L)
//        {
//            j = j - L;
//            L = L / 2;
//        }
//        j = j + L;
//        if (i<j)
//        {
//            tempReal = outputReal[i];
//            tempImaginary = outputImag[i];
//            outputReal[i] = outputReal[j];
//            outputImag[i] = outputImag[j];
//            outputReal[j] = tempReal;
//            outputImag[j] = tempImaginary;
//        }
//    }
//    L = 0;
//    m = 1;
//    n = k / 2;
//    //computation
//    for (i = k; i>1; i = (i >> 1))
//    {
//        L = m;
//        m = 2 * m;
//        o = 0;
//        for (j = 0; j<L; j++)
//        {
//            cos = cosine[o];
//            sin = sine[o];
//            o = o + n;
//            for (p = j; p<k; p = p + m)
//            {
//                q = p + L;
//                xt = cos*outputReal[q] - sin*outputImag[q];
//                yt = sin*outputReal[q] + cos*outputImag[q];
//                outputReal[q] = (outputReal[p] - xt);
//                outputImag[q] = (outputImag[p] - yt);
//                outputReal[p] = (outputReal[p] + xt);
//                outputImag[p] = (outputImag[p] + yt);
//            }
//        }
//        n = n >> 1;
//    }
//}
//
//void overlapAdd::IFFT(const float* inputReal, const float* inputImag, float *outputReal, float *outputImag, const int nFFT)
//{
//    float arg;
//    for (int i = 0; i<nFFT / 2; i++)
//    {
//        arg = -2 * M_PI*i / nFFT;
//        cosine[i] = cos(arg);
//        sine[i] = sin(arg);
//    }
//    int i, j, k, L, m, n, o, p, q;
//    float tempReal, tempImaginary, cos, sin, xt, yt;
//    k = nFFT;
//    for (i = 0; i<k; i++)
//    {
//        outputReal[i] = inputReal[i];
//        outputImag[i] = (-1)*inputImag[i];
//    }
//
//    j = 0;
//    m = k / 2;
//    //bit reversal
//    for (i = 1; i<(k - 1); i++)
//    {
//        L = m;
//        while (j >= L)
//        {
//            j = j - L;
//            L = L / 2;
//        }
//        j = j + L;
//        if (i<j)
//        {
//            tempReal = outputReal[i];
//            tempImaginary = outputImag[i];
//            outputReal[i] = outputReal[j];
//            outputImag[i] = outputImag[j];
//            outputReal[j] = tempReal;
//            outputImag[j] = tempImaginary;
//        }
//    }
//    L = 0;
//    m = 1;
//    n = k / 2;
//    //computation
//    for (i = k; i>1; i = (i >> 1))
//    {
//        L = m;
//        m = 2 * m;
//        o = 0;
//        for (j = 0; j<L; j++)
//        {
//            cos = cosine[o];
//            sin = sine[o];
//            o = o + n;
//            for (p = j; p<k; p = p + m)
//            {
//                q = p + L;
//                xt = cos*outputReal[q] - sin*outputImag[q];
//                yt = sin*outputReal[q] + cos*outputImag[q];
//                outputReal[q] = (outputReal[p] - xt);
//                outputImag[q] = (outputImag[p] - yt);
//                outputReal[p] = (outputReal[p] + xt);
//                outputImag[p] = (outputImag[p] + yt);
//            }
//        }
//        n = n >> 1;
//    }
//    for (i = 0; i<k; i++)
//    {
//        outputReal[i] = outputReal[i] / k;
//        outputImag[i] = outputImag[i] / k;
//    }
//}
