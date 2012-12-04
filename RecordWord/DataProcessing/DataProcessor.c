//
//  DataProcessor.c
//  AudioRecorder
//
//  Created by Sammy Long on 12/9/28.
//  Copyright (c) 2012年 Apexlearn Inc. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#import <Accelerate/Accelerate.h>

#import "alDCT.h"

static const float cPreEmCoef = 0.97;

#define	N   256
#define LOG2N 8
#define N_2 128

static const size_t cBitRate = 16;
static const int cSampleRate = 8000;
static const size_t cFrameShift = 80;
static const size_t cNumChans = 26;
static const size_t cNumCeps = 12;
static const size_t cCepLifter = 22;
static const uint8_t d_Ascii = 100;
static const uint8_t a_Ascii = 97;
static const uint8_t t_Ascii = 116;
static const uint8_t cMFCCDim = 39;// include E, D, A

void moveToDataChunkAndNormalize(uint8_t *wave, size_t fileSize, float **data, size_t *size) {
    
    // find 'data' section
    int i=0;
    for (i=0; i<fileSize; i++) {
        if (wave[i] == d_Ascii) {
            if (wave[i+1] == a_Ascii && wave[i+2] == t_Ascii && wave[i+3] == a_Ascii) {
                printf("hooray\n");
                break;
            }
        }
    }
    // put size
    size_t *sizePtr = (size_t *)(wave + i + 4);
    *size = (*sizePtr)/(cBitRate/8);
    printf("data chunk size %zd\n", *size);
    size_t numberOfSamples = *size;
    
    float numFramesf = (float)(numberOfSamples - N)/(float)cFrameShift;
    int numFrames = (int)ceil(numFramesf);
    size_t numSampPaddedZero = numFrames * cFrameShift + N;
    // update size
    *size = numSampPaddedZero;
    // put samples
    int16_t *dataPtr = (int16_t *)(wave + i + 4 + 4);
    float* samples = malloc(sizeof(float_t)*numSampPaddedZero);
    for (i=0; i<numberOfSamples; i++) {
        samples[i] = (float)dataPtr[i]/32767.0;
    }
    // padded zero
    for (i=numberOfSamples; i<numSampPaddedZero; i++) {
        samples[i] = 0.0f;
    }
    *data = samples;

}

float mel(int k,float fres)
{
    return 1127 * log(1 + k*fres);
}

void initFBank(float ** cfPtr, size_t *numChansPlus2, int8_t ** loChanPtr, float **loWtPtr ) {
    float *cf = (float *)malloc(sizeof(float)*(cNumChans+2));
    if (cf == NULL) {
        printf("Error malloc cf.\n ");
        return;
    }
    *cfPtr = cf;
    *numChansPlus2 = cNumChans + 2;
    int maxChan = cNumChans + 1;
    float mlo, mhi, ms;
    mlo = 0.0;
    float fres = (float)cSampleRate/(float)N;
    fres/= 700.0f;
    mhi = mel(N_2, fres);
    ms = mhi - mlo;
    cf[0] = 0.0f;
    for (int chan=1; chan <= maxChan; chan++) {
        cf[chan] = ((float)chan/(float)maxChan)*ms + mlo;
    }
    int8_t *loChan = (int8_t *)malloc(sizeof(int8_t)*N_2);
    if (loChan == NULL) {
        printf("Error malloc loChan.\n ");
        return;
    }
    *loChanPtr = loChan;
    int klo, khi;
    klo = 1; khi = N_2 - 1;
    for (int k=0,chan=1; k<N_2; k++){
        float melk = mel(k,fres);
        if (k<klo || k>khi) loChan[k]=-1;
        else {
            while (cf[chan] < melk  && chan<=maxChan) ++chan;
            loChan[k] = chan-1;
        }
    }
    float *loWt = (float *)malloc(sizeof(float_t)*N_2);
    if (loWt == NULL) {
        printf("Error malloc loWt.\n ");
        return;
    }
    *loWtPtr = loWt;

    for (int k=0; k<N_2; k++) {
        int chan = loChan[k];
        if (k<klo) loWt[k]=0.0;
        else {
            if (chan>0) {
                loWt[k] = ((cf[chan+1] - mel(k,fres)) /
                              (cf[chan+1] - cf[chan]));
            } else {
                loWt[k] = (cf[1]-mel(k,fres))/(cf[1] - mlo);
            }
        }
    }    
}

void preEmphasis(float *s, size_t size) {
    for (int i= size -1 ;i>=1;i--) {
        s[i] -= s[i-1]*cPreEmCoef;
    }
    s[0] *= 1.0-cPreEmCoef;
}

void wave2fbank(float *segment, float *m, int8_t *loChan, float *loWt) {
    // clear channel
    for (int i=0; i<cNumChans; i++) {
        m[i] = 0.0f;
    }
    // fill bins
    for (int k = 1; k<N_2 ; k++) {
        int bin = loChan[k];
        float t1 = loWt[k]*segment[k];
        if (bin>0) {
            m[bin - 1] += t1;
        }
        if (bin<cNumChans) {
            m[bin] += segment[k] - t1;
        }
    }
    for (int bin=0; bin<cNumChans; bin++) {
        m[bin] = logf(m[bin]);
    }
}

void fbank2mfcc(float *m, float *c) {
    float pi_factor = M_PI/(float)cNumChans;
    float mfnorm = sqrt(2.0/(float)cNumChans);
    for (int j=0; j<cNumCeps; j++)  {
        c[j] = 0.0;
        float x = (float)(j+1) * pi_factor;
        for (int k=0; k<cNumChans; k++) {
            c[j] += m[k] * cos(x*((k+1)-0.5));
        }
        c[j] *= mfnorm;
    }
}

void genCepWin(float *cepWin) {
    float a, Lby2;
    a = M_PI/cCepLifter;
    Lby2 = (float)cCepLifter/2.0;
    for (int i=0;i<cNumCeps;i++) {
        cepWin[i] = 1.0 + Lby2*sin((i+1) * a);
    }
}

void weightCepstrum (float * c, float *cepWin) {
    for (int i=0;i<cNumCeps;i++) {
        c[i] *= cepWin[i];
    }
}
/*
void regress (float *c, float *d) {
    static const int delwin = 2;
    float sigmaT2 = 0.0f;
    for (int t=1;t<=delwin;t++) {
        sigmaT2 += t*t;
    }
    sigmaT2 *= 2.0f;
    for (int i=2;i<cNumCeps-2;i++) {
        float sum = 0.0f;
        for (int t=1; t<=delwin; t++) {
            sum += t*(c[i+t] - c[i-t]);
        }
        d[i] = sum/sigmaT2;
    }
}
*/
// feature_extraction
void MFCC(uint8_t *wave, size_t wave_size) {
    printf("wave size %zd\n", wave_size);
    float *samples = NULL;
    size_t size = 0;
    
    // normalization
    moveToDataChunkAndNormalize(wave, wave_size, &samples, &size);
    
    // create ham window lookup table
    float *hamWin = (float *)malloc(N * sizeof(float_t));
    if (hamWin == NULL) {
        printf("Error malloc hamming window lookup table.\n");
        return;
    }
    vDSP_hamm_window(hamWin, N, 0);
    
    // initFBank
    float *cf = NULL;
    size_t cfSize = 0;
    int8_t *loChan = NULL;
    float *loWt = NULL;
    initFBank(&cf, &cfSize, &loChan, &loWt);
    
    // create fft setup
    FFTSetup setup = vDSP_create_fftsetup(LOG2N, FFT_RADIX2);
    DSPSplitComplex fftData;
    fftData.realp = (float *) malloc(N_2 * sizeof(float_t));
    fftData.imagp = (float *) malloc(N_2 * sizeof(float_t));
    if (fftData.realp == NULL || fftData.imagp == NULL) {
        printf("Error malloc fftData.\n ");
        return;
    }
    // create buffer for each frame
    float *segment = malloc(sizeof(float) * N);
    if (segment == NULL) {
        printf("Error malloc segment.\n ");
        return;
    }
    // buffer for m (filterbank channels)
    float *m = (float *)malloc(sizeof(float_t) * cNumChans);
    if (m == NULL) {
        printf("Error malloc m.\n ");
        return;
    }
    int numFrames = (size - N)/cFrameShift;

    // buffer for all mfccs
    float *mfcc = (float *)malloc(sizeof(float_t) * cMFCCDim * (numFrames));
    if (mfcc == NULL) {
        printf("Error malloc c.\n ");
        return;
    }
    // create cepWin
    float *cepWin = (float *)malloc(sizeof(float_t) * cNumCeps);
    if (cepWin == NULL) {
        printf("Error malloc cepWin.\n ");
        return;
    }
    genCepWin(cepWin);


    // grab the first frame (for now)
    for (int i=0; i<numFrames; i++) {
        memcpy(segment, samples + cFrameShift*i, sizeof(float_t) * N);
        // PreEmphasis
        preEmphasis(segment, N);
        // ham
        vDSP_vmul(segment, 1, hamWin, 1, segment, 1, N);
        // fft
        // split input with even indices to real and odd indices to imag
        vDSP_ctoz((DSPComplex *)segment, 2, &fftData, 1, N_2);
        // in-place real fft 
        vDSP_fft_zrip(setup, &fftData, 1, LOG2N, FFT_FORWARD);
/*
        for (int j=0; j<N_2; j++) {
            printf("%d: real %f, imag %f\n", j, fftData.realp[j], fftData.imagp[j]);
        }
 */
        // get mag
        for(int j=0; j<N_2; j++) {
            segment[j] = sqrt(pow(fftData.realp[j], 2) + pow(fftData.imagp[j], 2));
            //printf("%d: %f\n", j, segment[j]);
        }
        wave2fbank(segment, m, loChan, loWt);
        fbank2mfcc(m, mfcc + i*cMFCCDim);
        // weight cepstrum
        weightCepstrum(mfcc + i*cMFCCDim, cepWin);
    }
    
    //regress
    // first derivative
    
    // denominator
    float denom = 0.0f;
    for (int t=1; t<=2; t++) {
        denom += t*t;
    }
    denom *= 2.0f;
    float *prev = NULL;
    float *next = NULL;
    for (int i=0; i<numFrames; i++) {
        for (int j=0; j<cNumCeps; j++) {
            float sum = 0.0f;
            next = prev = mfcc + i * cMFCCDim + j;
            for (int t=1; t<=2; t++) {
                if (i+t<numFrames) next += cMFCCDim;
                if (i-t>=0) prev -= cMFCCDim;
                sum += t * ( *next - *prev);
            }
            sum /= denom;
            mfcc[i * cMFCCDim + j + cNumCeps + 1] = sum;
        }
    }

    //debug print
    for (int i=191; i<numFrames; i++) {
        printf("frame %d\n", i);
        for (int j=0; j<26; j++) {
            printf("%f, ", mfcc[i * cMFCCDim + j]);
        }
        printf("\n");
    }
    
    // free allocated memory
    vDSP_destroy_fftsetup(setup);
    free(fftData.realp);
    free(fftData.imagp);
    free(segment);
    free(hamWin);
    free(cf);
    free(loChan);
    free(loWt);
    free(m);
    free(mfcc);
    free(samples);
}

void FourierTransformTest(void) {
    float * OriginalReal = (float *)malloc(N * sizeof(float));
    for (int i = 0; i < N ; i++) {
        OriginalReal[i] = (float) (i + 1);
    }
    // create fft setup
    FFTSetup setup = vDSP_create_fftsetup(LOG2N, FFT_RADIX2);
    DSPSplitComplex fftData;
    fftData.realp = (float *) malloc(N_2 * sizeof(float));
    fftData.imagp = (float *) malloc(N_2 * sizeof(float));
    // split input with even indices to real and odd indices to imag
    vDSP_ctoz((DSPComplex *)OriginalReal, 2, &fftData, 1, N_2);
/*
    printf("after ctoz:");
    for (int i = 0; i < N_2; i++) {
        printf("i=%d: %.2f + %.2fi\n",i ,fftData.realp[i], fftData.imagp[i]);
    }
*/
    // FFT real in-place
    vDSP_fft_zrip(setup, &fftData, 1, LOG2N, FFT_FORWARD);

    printf("after zrip forward:");
    for (int i = 0; i < N_2; i++) {
        printf("i=%d: %.2f + %.2fi\n",i ,fftData.realp[i], fftData.imagp[i]);
    }

    vDSP_fft_zrip(setup, &fftData, 1, LOG2N, FFT_INVERSE);
    float scale = 1.0f/(2.0f * N);
    vDSP_vsmul(fftData.realp, 1, &scale, fftData.realp, 1, N_2);
    vDSP_vsmul(fftData.imagp, 1, &scale, fftData.imagp, 1, N_2);
/*
    printf("after zrip inverse:");
    for (int i = 0; i < N_2; i++) {
        printf("i=%d: %.2f + %.2fi\n",i ,fftData.realp[i], fftData.imagp[i]);
    }
*/
    // destroy fft setup
    vDSP_destroy_fftsetup(setup);
    // free allocated memory
    free(OriginalReal);
    free(fftData.realp);
    free(fftData.imagp);
}
