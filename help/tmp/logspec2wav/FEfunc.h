/*---------------------------------------------------------------------------
 *
 * Mel Cepstrum Parametrisation (Track 1 Front End) v2.0 for Distributed
 * Speech Recognition (DSR), basic front-end functions
 *
 * FILE NAME: FEfunc.h
 * PURPOSE:   This function package contains basic front-end functions. All
 *            functions except the initializations are called frame by frame
 *            from main
 *
 * This software has been released to the members of the Aurora Project
 * to be used in the validation and verification of the proposed ETSI
 * DSR Track 1 Front End.
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/

#ifndef _FEFUNC_H
#define _FEFUNC_H

/*----------------------
 * Constant Definitions
 *----------------------*/
#define M_PI        3.14159265358979323846
#define M_SQRT2     1.41421356237309504880

#define PI          M_PI
#define PIx2        6.28318530717958647692

/*-----------------
 * Data Structures
 *-----------------*/
typedef struct FFT_Window_tag
{
    int StartingPoint;
    int Length;
    float *Data;
    struct FFT_Window_tag *Next;
}
FFT_Window;                     /* Structure for FFT window (one triangle in
				   the chained list) */

/*---------------------
 * Function Prototypes
 *---------------------*/
void DCOffsetFilter( float *CircBuff, int BSize, int *BPointer, int nSamples );
void InitializeHamming (float *win, int len);
void Window (float *data, float *win, int len);
void DeWindow (float *data, float *win, int len);
void rfft (float *x, int n, int m);
void rifft (float *x, int n, int m);
void InitFFTWindows (FFT_Window * FirstWin,
                     float StFreq,
                     float SmplFreq,
                     int FFTLength,
                     int NumChannels);
void ReleaseFFTWindows (FFT_Window *FirstWin );
void ComputeTriangle (FFT_Window * FirstWin);
void MelFilterBank (float *SigFFT, FFT_Window * FirstWin);

float *InitDCTMatrix (int NumCepstralCoeff, int NumChannels);
void DCT (float *Data, float *Mx, int NumCepstralCoeff, int NumChannels);

#endif
