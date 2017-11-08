/*---------------------------------------------------------------------------
 *
 * Mel Cepstrum Parametrisation (Track 1 Front End) v2.0 for Distributed
 * Speech Recognition (DSR), basic front-end functions
 *
 * FILE NAME: FEfunc.c
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

#define _CRT_SECURE_NO_DEPRECATE

/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdlib.h>
#include <math.h>
#include "FEfunc.h"


/*---------------------------------------------------------------------------
 * FUNCTION NAME: DCOffsetFilter
 *
 * PURPOSE:       DC offset removal from speech waveform
 *
 * INPUT:
 *   CircBuff     Pointer to input circular buffer
 *   BSize        Buffer size
 *   BPointer     Pointer to buffer pointer
 *   nSamples     Number of samples to be filtered
 *
 * OUTPUT
 *                Filtered data in the circular buffer
 *                last output sample pointed by the buffer pointer
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void
DCOffsetFilter( float *CircBuff, int BSize, int *BPointer, int nSamples )
{
  int i;

  for ( i=0; i<nSamples; i++ )
    {
      /* y[n]=x[n]-x[n-1]+0.999*y[n-1] */
      CircBuff[(*BPointer+1)%BSize] = (float) (CircBuff[(*BPointer+2)%BSize] - CircBuff[(*BPointer+1)%BSize]+0.999*CircBuff[*BPointer]);
      *BPointer=(*BPointer+1)%BSize;
    }
}


/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitializeHamming
 *
 * PURPOSE:       Initializes Hamming window coefficients
 *
 * INPUT:
 *   win          Pointer to window buffer
 *   len          Window length
 *
 * OUTPUT
 *                Hamming window coefficients stored in window buffer pointed
 *                to by *win*
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
InitializeHamming (float *win, int len)
{
    int i;

    for (i = 0; i < len / 2; i++)
        win[i] = (float) (0.54 - 0.46 * cos (PIx2 * i / (len - 1)));
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: Window
 *
 * PURPOSE:       Performs windowing on input speech frame (multiplies input
 *                samples by the corresponding window coefficients)
 *
 * INPUT:
 *   data         Pointer to input speech buffer
 *   win          Pointer to window buffer
 *   len          Window (or frame) length
 *
 * OUTPUT
 *                Windowed speech frame stored at the same place as the
 *                original speech samples (pointed by *data*)
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
Window (float *data, float *win, int len)
{
    long i;

    for (i = 0; i < len / 2; i++)
        data[i] *= win[i];
    for (i = len / 2; i < len; i++)
        data[i] *= win[len - 1 - i];
}

void 
DeWindow (float *data, float *win, int len)
{
	long i;

	for (i = 0; i < len / 2; i++)
		data[i] /= win[i];
	for (i = len / 2; i < len; i++)
		data[i] /= win[len - 1 - i];
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: rfft
 *
 * PURPOSE:       Real valued, in-place split-radix FFT
 *
 * INPUT:
 *   x            Pointer to input and output array
 *   n            Length of FFT, must be power of 2
 *
 * OUTPUT         Output order
 *                  Re(0), Re(1), ..., Re(n/2), Im(N/2-1), ..., Im(1)
 *
 * RETURN VALUE
 *   none
 *
 * DESIGN REFERENCE:
 *                IEEE Transactions on Acoustic, Speech, and Signal Processing,
 *                Vol. ASSP-35. No. 6, June 1987, pp. 849-863.
 *
 *                Subroutine adapted from fortran routine pp. 858-859.
 *                Note corrected printing errors on page 859:
 *                    SS1 = SIN(A3) -> should be SS1 = SIN(A);
 *                    CC3 = COS(3)  -> should be CC3 = COS(A3)
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
rfft (float *x, int n, int m)
{
    int j, i, k, is, id;
    int i0, i1, i2, i3, i4, i5, i6, i7, i8;
    int n2, n4, n8;
    float xt, a0, e, a, a3;
    float t1, t2, t3, t4, t5, t6;
    float cc1, ss1, cc3, ss3;
    float *r0;

    /* Digit reverse counter */

    j = 0;
    r0 = x;

    for (i = 0; i < n - 1; i++)
    {

        if (i < j)
        {
            xt = x[j];
            x[j] = *r0;
            *r0 = xt;
        }
        r0++;

        k = n >> 1;

        while (k <= j)
        {
            j = j - k;
            k >>= 1;
        }
        j += k;
    }

    /* Length two butterflies */
    is = 0;
    id = 4;

    while (is < n - 1)
    {

        for (i0 = is; i0 < n; i0 += id)
        {
            i1 = i0 + 1;
            a0 = x[i0];
            x[i0] += x[i1];
            x[i1] = a0 - x[i1];
        }

        is = (id << 1) - 2;
        id <<= 2;
    }

    /* L shaped butterflies */
    n2 = 2;
    for (k = 1; k < m; k++)
    {
        n2 <<= 1;
        n4 = n2 >> 2;
        n8 = n2 >> 3;
        e = (float) ((M_PI * 2) / n2);
        is = 0;
        id = n2 << 1;
        while (is < n)
        {
            for (i = is; i <= n - 1; i += id)
            {
                i1 = i;
                i2 = i1 + n4;
                i3 = i2 + n4;
                i4 = i3 + n4;
                t1 = x[i4] + x[i3];
                x[i4] -= x[i3];
                x[i3] = x[i1] - t1;
                x[i1] += t1;

                if (n4 != 1)
                {
                    i1 += n8;
                    i2 += n8;
                    i3 += n8;
                    i4 += n8;
                    t1 = (float) ((x[i3] + x[i4]) / M_SQRT2);
                    t2 = (float) ((x[i3] - x[i4]) / M_SQRT2);
                    x[i4] = x[i2] - t1;
                    x[i3] = -x[i2] - t1;
                    x[i2] = x[i1] - t2;
                    x[i1] = x[i1] + t2;
                }
            }
            is = (id << 1) - n2;
            id <<= 2;
        }

        for (j = 1; j < n8; j++)
        {
            a = j * e;
            a3 = 3 * a;
            cc1 = (float) cos (a);
            ss1 = (float) sin (a);
            cc3 = (float) cos (a3);
            ss3 = (float) sin (a3);

            is = 0;
            id = n2 << 1;

            while (is < n)
            {
                for (i = is; i <= n - 1; i += id)
                {
                    i1 = i + j;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    i4 = i3 + n4;
                    i5 = i + n4 - j;
                    i6 = i5 + n4;
                    i7 = i6 + n4;
                    i8 = i7 + n4;
                    t1 = x[i3] * cc1 + x[i7] * ss1;
                    t2 = x[i7] * cc1 - x[i3] * ss1;
                    t3 = x[i4] * cc3 + x[i8] * ss3;
                    t4 = x[i8] * cc3 - x[i4] * ss3;
                    t5 = t1 + t3;
                    t6 = t2 + t4;
                    t3 = t1 - t3;
                    t4 = t2 - t4;
                    t2 = x[i6] + t6;
                    x[i3] = t6 - x[i6];
                    x[i8] = t2;
                    t2 = x[i2] - t3;
                    x[i7] = -x[i2] - t3;
                    x[i4] = t2;
                    t1 = x[i1] + t5;
                    x[i6] = x[i1] - t5;
                    x[i1] = t1;
                    t1 = x[i5] + t4;
                    x[i5] = x[i5] - t4;
                    x[i2] = t1;
                }
                is = (id << 1) - n2;
                id <<= 2;
            }
        }
    }
}


void 
rifft (float *x, int n, int m)
{
    int j, i, k, is, id;
    int i0, i1, i2, i3, i4, i5, i6, i7, i8;
    int n2, n4, n8;
    float xt, a0, e, a, a3;
    float t1, t2, t3, t4, t5;
    float cc1, ss1, cc3, ss3;
    float *r0;

    /* L shaped butterflies */
    n2 = 2*n;
    for (k = 0; k < m-1; k++)
    {
		is = 0; 
		id = n2;
        n2 >>= 1;
        n4 = n2 >> 2;
        n8 = n4 >> 1;
        e = (float) ((M_PI * 2) / n2);
        
        while (is < n-1)
        {
            for (i = is; i <= n - 1; i += id)
            {
                i1 = i;
                i2 = i1 + n4;
                i3 = i2 + n4;
                i4 = i3 + n4;
                t1 = x[i1] - x[i3];
                x[i1] += x[i3];
				x[i2] = x[i2]*2;
                x[i3] = t1 - 2*x[i4];
                x[i4] = t1 + 2*x[i4];

                if (n4 != 1)
                {
                    i1 += n8;
                    i2 += n8;
                    i3 += n8;
                    i4 += n8;
                    t1 = (float) ((x[i2] - x[i1]) / M_SQRT2);
                    t2 = (float) ((x[i4] + x[i3]) / M_SQRT2);
                    x[i1] += x[i2];
                    x[i2] = x[i4] - x[i3];
                    x[i3] = -2*(t1 + t2);
                    x[i4] = 2*(t1 - t2);
                }
            }
            is = (id << 1) - n2;
            id <<= 2;
        }

        for (j = 1; j < n8; j++)
        {
            a = j * e;
            a3 = 3 * a;
            cc1 = (float) cos (a);
            ss1 = (float) sin (a);
            cc3 = (float) cos (a3);
            ss3 = (float) sin (a3);

            is = 0;
            id = n2 << 1;

            while (is < n-1)
            {
                for (i = is; i <= n - 1; i += id)
                {
                    i1 = i + j;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    i4 = i3 + n4;
                    i5 = i + n4 - j;
                    i6 = i5 + n4;
                    i7 = i6 + n4;
                    i8 = i7 + n4;
					t1 = x[i1] - x[i6];
					x[i1] += x[i6];
					t2 = x[i5] -x[i2];
					x[i5] += x[i2];
					t3 = x[i8] + x[i3];
					x[i6] = x[i8] - x[i3];
					t4 = x[i4] + x[i7];
					x[i2] = x[i4] - x[i7];
					t5 = t1 - t4;
					t1 += t4;
					t4 = t2 - t3;
					t2 += t3;
					x[i3] = t5*cc1 + t4*ss1;
					x[i7] = -t4*cc1 + t5*ss1;
					x[i4] = t1*cc3 - t2*ss3;
					x[i8] = t2*cc3 + t1*ss3;
                }
                is = (id << 1) - n2;
                id <<= 2;
            }
        }
    }

	 /* Length two butterflies */
    is = 0;
    id = 4;

    while (is < n - 1)
    {

        for (i0 = is; i0 < n; i0 += id)
        {
            i1 = i0 + 1;
            a0 = x[i0];
            x[i0] += x[i1];
            x[i1] = a0 - x[i1];
        }

        is = (id << 1) - 2;
        id <<= 2;
    }

	/* Digit reverse counter */

    j = 0;
    r0 = x;

    for (i = 0; i < n - 1; i++)
    {

        if (i < j)
        {
            xt = x[j];
            x[j] = *r0;
            *r0 = xt;
        }
        r0++;

        k = n >> 1;

        while (k <= j)
        {
            j = j - k;
            k >>= 1;
        }
        j += k;
    }

	for (i = 0; i < n; i++)
    {
		x[i] /= n;
	}

}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitFFTWindows
 *
 * PURPOSE:       Initializes data structure for FFT windows (mel filter bank).
 *                Computes starting point and length of each window, allocates
 *                memory for window coefficients.
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *   StFreq       Starting frequency of mel filter bank
 *   SmplFreq     Sampling frequency
 *   FFTLength    FFT length
 *   NumChannels  Number of channels
 *
 * OUTPUT
 *                Chained list of FFT window data structures. NOTE FFT window
 *                coefficients are not computed yet.
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
InitFFTWindows (FFT_Window * FirstWin,
                float StFreq,
                float SmplFreq,
                int FFTLength,
                int NumChannels)

{
    int i, TmpInt;
    float freq, start_mel, fs_per_2_mel;
    FFT_Window *p1, *p2;

    /* Constants for calculation */
    start_mel = (float) (2595.0 * log10 (1.0 + (float) StFreq / 700.0));
    fs_per_2_mel = (float) (2595.0 * log10 (1.0 + (SmplFreq / 2) / 700.0));

    p1 = FirstWin;

    for (i = 0; i < NumChannels; i++)
    {
        /* Calculating mel-scaled frequency and the corresponding FFT-bin */
        /* number for the lower edge of the band                          */
        freq = (float) (700 * (pow (10, (start_mel + (float) i / (NumChannels + 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0));
        TmpInt = (int) (FFTLength * freq / SmplFreq + 0.5);

        /* Storing */
        p1->StartingPoint = TmpInt;

        /* Calculating mel-scaled frequency for the upper edge of the band */
        freq = (float) (700 * (pow (10, (start_mel + (float) (i + 2) / (NumChannels + 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0));

        /* Calculating and storing the length of the band in terms of FFT-bins*/
        p1->Length = (int) (FFTLength * freq / SmplFreq + 0.5) - TmpInt + 1;

        /* Allocating memory for the data field */
        p1->Data = (float *) malloc (sizeof (float) * p1->Length);

        /* Continuing with the next data structure or close the last structure
	   with NULL */
        if (i < NumChannels - 1)
        {
            p2 = (FFT_Window *) malloc (sizeof (FFT_Window));
            p1->Next = p2;
            p1 = p2;
        }
        else
            p1->Next = NULL;
    }
    return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: ReleaseFFTWindows
 *
 * PURPOSE:       Releases memory allocated for FFT windows
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/

void ReleaseFFTWindows (FFT_Window *FirstWin )
{
  FFT_Window *p;

  while ( FirstWin->Next!=NULL )
    {
    p=FirstWin->Next->Next;
    free(FirstWin->Next->Data);
    free(FirstWin->Next);
    FirstWin->Next=p;
    }
  free(FirstWin->Data);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: ComputeTriangle
 *
 * PURPOSE:       Computes and stores FFT window coefficients (triangle points)
 *                into initialized chained list of FFT window structures
 *
 * INPUT:
 *   FirstWin     Pointer to first FFT window structure
 *
 * OUTPUT
 *                Chained list of FFT window data structures with correct
 *                window coefficients
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
ComputeTriangle (FFT_Window * FirstWin)
{
    FFT_Window *p1;

    int low_part_length, hgh_part_length, TmpInt=0, i, j;

    p1 = FirstWin;
    j = 0;
    while (p1)
    {
        low_part_length = p1->Next ?
	                  p1->Next->StartingPoint - p1->StartingPoint + 1 :
                          TmpInt - p1->StartingPoint + 1;
        hgh_part_length = p1->Length - low_part_length + 1;

        /* Lower frequency part of the triangle */
        for (i = 0; i < low_part_length; i++)
            p1->Data[i] = (float) (i + 1) / low_part_length;

        /* Higher frequency part of the triangle */
        for (i = 1; i < hgh_part_length; i++)
            p1->Data[low_part_length + i - 1] = (float) (hgh_part_length - i) /
	      hgh_part_length;

        /* Store upper edge (for calculating the last triangle) */
        TmpInt = p1->StartingPoint + p1->Length - 1;

        /* Next triangle ... */
        p1 = p1->Next;
    }
    return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: MelFilterBank
 *
 * PURPOSE:       Performs mel filtering on FFT magnitude spectrum using the
 *                filter bank defined by a chained list of FFT window
 *                structures
 *
 * INPUT:
 *   SigFFT       Pointer to signal FFT magnitude spectrum
 *   FirstWin     Pointer to the first channel of the filter bank (first
 *                element in the chained list of FFT window data structures)
 *
 * OUTPUT
 *                Filter bank outputs stored at the beginning of input signal
 *                FFT buffer pointed by *SigFFT*
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
MelFilterBank (float *SigFFT, FFT_Window * FirstWin)
{
    FFT_Window *p1;
    float Sum;
    int i, j;

    p1 = FirstWin;
    j = 0;
    while (p1)
    {
        Sum = 0.0;
        for (i = 0; i < p1->Length; i++)
            Sum += SigFFT[p1->StartingPoint + i] * p1->Data[i];
        SigFFT[j] = Sum;
        j++;
        p1 = p1->Next;
    }
    return;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: InitDCTMatrix
 *
 * PURPOSE:       Initializes matrix for DCT computation (DCT is implemented
 *                as matrix-vector multiplication). The DCT matrix is of size
 *                (NumCepstralCoeff-1)-by-NumChannels. The zeroth cepstral
 *                coefficient is computed separately (needing NumChannels
 *                additions and only one multiplication), so the zeroth row
 *                of DCT matrix corresponds to the first DCT basis vector, the
 *                first one to the second one, and so on up to
 *                NumCepstralCoeff-1.
 *
 * INPUT:
 *   NumCepstralCoeff
 *                Number of cepstral coeffficients
 *   NumChannels  Number of filter bank channels
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *                Pointer to the initialized DCT matrix
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
float *
InitDCTMatrix (int NumCepstralCoeff, int NumChannels)
{
    int i, j;
    float *Mx;

    /* Allocating memory for DCT-matrix */
    Mx = (float *) malloc (sizeof (float) * (NumCepstralCoeff - 1) *
			   NumChannels);

    /* Computing matrix entries */
    for (i = 1; i < NumCepstralCoeff; i++)
        for (j = 0; j < NumChannels; j++)
            Mx[(i - 1) * NumChannels + j] = (float) cos (PI * (float) i / (float) NumChannels * ((float) j + 0.5));
    return Mx;
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: DCT
 *
 * PURPOSE:       Computes DCT transformation of filter bank outputs, results
 *                in cepstral coefficients. The DCT transformation is
 *                implemented as matrix-vector multiplication. The zeroth
 *                cepstral coefficient is computed separately and appended.
 *                Final cepstral coefficient order is c1, c2, ...,c12, c0. The
 *                output is stored right after the input values in the memory.
 *                Since the mel filter bank outputs are stored at the beginning
 *                of the FFT magnitude array it shouldn`t cause any problems.
 *                Some memory saving can be done this way.
 *
 * INPUT:
 *   Data         Pointer to input data buffer (filter bank outputs)
 *   Mx           Pointer to DCT matrix
 *   NumCepstralCoeff
 *                Number of cepstral coefficients
 *   NumChannels  Number of filter bank channels
 *
 * OUTPUT
 *                Cepstral coefficients stored after the input filter bank
 *                values pointed to by *Data*
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
DCT (float *Data, float *Mx, int NumCepstralCoeff, int NumChannels)
{
    int i, j;

    /* Computing c1..c/NumCepstralCoeff-1/, storing result after the incoming
       data vector */
    for (i = 1; i < NumCepstralCoeff; i++)
    {
        Data[NumChannels + (i - 1)] = 0.0;
        for (j = 0; j < NumChannels; j++)
            Data[NumChannels + (i - 1)] += Data[j]
	      * Mx[(i - 1) * NumChannels + j];
    }

    /* Computing c0, as the last element of output vector */
    Data[NumChannels + NumCepstralCoeff - 1] = 0.0;
    for (i = 0; i < NumChannels; i++)
        Data[NumChannels + NumCepstralCoeff - 1] += Data[i];
    return;
}
