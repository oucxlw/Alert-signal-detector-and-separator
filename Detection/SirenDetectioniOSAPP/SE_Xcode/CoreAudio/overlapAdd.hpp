//
//  overlapAdd.hpp
//  frameOverlap
//
//  Created by Shankar, Nikhil on 8/28/19.
//  Copyright © 2019 default. All rights reserved.
//

#ifndef overlapAdd_hpp
#define overlapAdd_hpp

#include <stdio.h>

class overlapAdd
{
public:
    void init();
    ~overlapAdd();
    void process(const float *input,float *output, const int frameCount);
    void FFT(const float* input, float *outputReal, float *outputImag, const int nFFT);
    void IFFT(const float* inputReal, const float* inputImag, float *outputReal, float *outputImag, const int nFFT);
    void siren_det(const float *input, float *class_out, const int frameCount);
    void setFrameSize(int frameSize){frameS = frameSize;}
    void setFftSize(int fftSize){nFFT = fftSize;}
    void setfs(int fs1){fs = fs1;}
private:
    float *in_prev;
    float *inputBuffer;
    int frameS;
    int nFFT;
    int fs;
    float *win, *invWin;
    float *xReal1, *xImag1, *yReal, *yImag;
    float *xReal2, *xImag2;
    float *cosine, *sine;
    float *outputOld;
    float *mag;
    float *noiseEstimatePrev;
    float *Pprev;
    float *meanP;
    float *meanP2;
    float alphaC;
    float alphaCPrev;
    float *alphaHat;
    float snr;
    float *P;
    float *beta;
    float *varP;
    float *invQeq;
    float *tilda_QD;
    float *tilda_QV;
    float *Bmin;
    float *Bminsub;
    float meanInv,  Bc;
    float *PBC;
    int *kMode;
    float *actmin;
    float *actminSub;
    int subwc;
    int *limflag;
    float *Pmin;
    float noiseSlopeMax;
    int *p;
    int *temp1, *temp2;
    float *noiseEstimate;
    float *noiseMean;
    float *noisePower;
    float *noiseMu;
    float *postSNR;
    float *aprioriSNR;
    float *prevOutput;
    float *logSigma;
    float vadDecision;
    int vad;
    int prevVad;
    int count1, count0;
    float *postSNR1;
    float *aprioriSNR1;
    float *gain;
    float *prevGain;
    float *ensig;
    float *finalSig;
    float *angle;
    float *de;
    float *dewin;
    // Siren Detect Parameters
    float *log_mag;
    
};
#endif /* overlapAdd_hpp */
