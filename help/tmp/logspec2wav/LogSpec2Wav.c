/*---------------------------------------------------------------------------
 *
 * Mel Cepstrum Parametrisation (Track 1 Front End) v2.0 for Distributed
 * Speech Recognition (DSR), command line parsing and main function
 *
 * FILE NAME: FrontEnd.c
 * PURPOSE:   This file contains the main part of DSR front-end. Speech samples
 *            are read from input waveform file frame by frame. Feature
 *            extraction is performed for each frame by calling basic
 *            functions of the basic FE function package (see FEfunc.c).
 *            Feature vectors are output to file in HTK format.
 *            Command line arguments are handled by a command line parsing
 *            function.
 *
 * This software has been released to the members of the Aurora Project
 * to be used in the validation and verification of the proposed ETSI
 * DSR Track 1 Front End.
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
 
 //YONG XU (yong.xu.ustc@gmail.com), 2015.12, modified, it can not be used for commercial use. If so, please contact me.
 
#define _CRT_SECURE_NO_DEPRECATE

/*-----------------
 * File Inclusions
 *-----------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fileio.h"
#include "FEfunc.h"

/*------------------------------------
 * Definition of Front-end Parameters 
 *------------------------------------*/
#define SAMPLING_FREQ_1           8  /*8kHz */
#define SAMPLING_FREQ_2          11  /*11kHz */
#define SAMPLING_FREQ_3          16  /*16kHz */

#define FRAME_LENGTH_1          256  /*32ms */
#define FRAME_LENGTH_2          256  /*23.27ms */
#define FRAME_LENGTH_3          512  /*32ms */

#define FRAME_SHIFT_1           128  /*16ms */
#define FRAME_SHIFT_2           110  /*10ms */
#define FRAME_SHIFT_3           256  /*16ms */

#define PRE_EMPHASIS           0.97

#define ENERGYFLOOR_FB        -50.0  /*0.0 */
#define ENERGYFLOOR_logE      -50.0  /*4.0 */

#define FFT_LENGTH_1            256
#define FFT_LENGTH_2            256
#define FFT_LENGTH_3            512
#define MAXWINDOWSIZE           512

#define NUM_CHANNELS             23
#define STARTING_FREQ_1        64.0  /*55.401825 */
#define STARTING_FREQ_2        64.0  /*63.817818 */
#define STARTING_FREQ_3        64.0  /*74.238716 */

#define NUM_CEP_COEFF            13  /* c1..c12 + c0 */

//#define MAX_VECTOR_NUM 			10000
//#define MAX_SAMPLE_NUM			100000
#define MAX_VECTOR_NUM 			200000 //2013.11.7因为要处理一个较长的实际语音，把它改大了点
#define MAX_SAMPLE_NUM			2000000

#define OLA_KIND				1 // 默认是1
#define POSTPROCESS				 1 // 默认是0
#define SMOOTHPROCESS			1 // 默认是0
//#define SMOOTH_WIN				1 // 默认的是1
#define SMOOTH_WIN				2 // 默认的是1
#define SWITCHPOINT			    36
#define THRESHOLD1				-2.1
#define THRESHOLD2				-3.43
#define NOISE_FRAME_NUM			10

#define LOWSEGSNR				-20
#define HIGHSEGSNR				30
#define DYNRANGE				50

/*-------------------------------
 * Global Definitions and Macros
 *-------------------------------*/
#define	BOOLEAN int
#define	FALSE 0
#define	TRUE (!FALSE)
#define WAIT_A_KEY while( !getchar() );
#define PRINTMOD 15
#define IEEE_LE 0
#define IEEE_BE 1

/*----------------------------------------------
 * Global Variable Definitions and Declarations
 *----------------------------------------------*/
BOOLEAN QuietMode = FALSE,      /* Suppress output to stderr */
  FsSpecified = FALSE,          /* Sampling frequency specified */
  SwapSpecified = FALSE,        /* Byte swap for raw data files specified */
  InputKindSpecified = FALSE,   /* Input file format specified */
  NoOutHeaderSpecified = FALSE, /* No output HTK header option specified */
  Noc0 = FALSE,                 /* No c0 coefficient to output feature vector */
  NologE = FALSE,               /* No logE component to output feature vector */ 
  NoisyInfo = FALSE;

FILE *fp_in_clean = NULL,       /* Input HTK, NIST or raw data file */
	 *fp_in_noisy = NULL,       /* Input HTK, NIST or raw data file */
     *fp_in_fea   = NULL,		/* Input enhanced feature file */
     *fp_out_info = NULL,		/* Output info file*/
	 *fp_out_noisy_info = NULL, /* Output noisy info file*/
     *fp_out_proc = NULL;       /* Output processed wave file*/

char InFilename[199],           /* Name of input file */
  InFilename_Fea[199],			/* Name of input feature file */
  InFilename_Noisy[199],		/* Name of input feature file */
  OutFilename_Info[199],		/* Name of output info file */
  OutFilename_Noisy_Info[199],  /* Name of output noisy info file */
  OutFilename_Proc[199],        /* Name of output processed wave file */
  InputKind[10] = "NIST";       /* Input file format */

int SamplingFrequency = 16000,  /* SamplingFrequency */
  NativeByteOrder = IEEE_LE,    /* Native byte ordering */
  InputByteOrder,               /* Default input byte ordering */
  OutParmKind=9;             /* Output parameter kind (MFCC, MFCC_E, MFCC_0, */
                                /* or MFCC_0_E as default) */

float out_proc_vectors[MAXWINDOWSIZE + 1][MAX_VECTOR_NUM];

float clean_spec_vectors[MAXWINDOWSIZE + 1][MAX_VECTOR_NUM];
float noisy_spec_vectors[MAXWINDOWSIZE + 1][MAX_VECTOR_NUM];
float denoise_spec_vectors[MAXWINDOWSIZE + 1][MAX_VECTOR_NUM];

int waveform_num;
float out_proc_wave[MAX_SAMPLE_NUM];
float waveform_count[MAX_SAMPLE_NUM];
short waveform[MAX_SAMPLE_NUM];

float segsnr;
float segsnr_noisy;
float lsd;
float lsd_noisy;

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ParseCommLine
 *
 * PURPOSE:       Parses command line arguments, opens input and output files
 *
 * INPUT:
 *   argc         Number of command line arguments
 *   argv         Array of command line arguments
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   FALSE        In case of any errors
 *   TRUE         Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
static BOOLEAN 
ParseCommLine (int argc, char *argv[])
{
    int mark = 0;

    while (argc)
    {
        if (strcmp (argv[mark], "-q") == 0)
        {
            QuietMode = TRUE;
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-fs") == 0)
        {
            FsSpecified = TRUE;
            --argc;
            ++mark;
            SamplingFrequency = 1000 * atoi (argv[mark]);
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-swap") == 0)
        {
            SwapSpecified = TRUE;
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-F") == 0)
        {
            InputKindSpecified = TRUE;
            --argc;
            ++mark;
            strcpy (InputKind, argv[mark]);
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-noh") == 0)
        {
            NoOutHeaderSpecified = TRUE;
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-noc0") == 0)
        {
            Noc0 = TRUE;
            --argc;
            ++mark;
        }
        else if (strcmp (argv[mark], "-nologE") == 0)
        {
            NologE = TRUE;
            --argc;
            ++mark;
        }
		else if (strcmp (argv[mark], "-ni") == 0)
		{
			NoisyInfo = TRUE;
			--argc;
			++mark;
		}
        else if (argv[mark][0] == '-')
        {
            fprintf (stderr, "WARNING:  Un-recognized flag '%s' !\r\n", argv[mark]);
            --argc;
            ++mark;
        }
        else
        {
            if (fp_out_proc)         /* Last argument string ERROR! */
            {
                fprintf (stderr, "ERROR:   Too many input arguments!\r\n");
                return FALSE;
            }
			else if(fp_out_info)  /* Fifth argument string (output processed wave file) */
			{
				strcpy (OutFilename_Proc, argv[mark]);
                fp_out_proc = fopen (OutFilename_Proc, "wb");
                if (fp_out_proc == NULL)
                {
                    fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", OutFilename_Proc);
                    return FALSE;
                }

			}
			else if(fp_in_fea)  /* Fourth argument string (output infomation file) */
			{
				strcpy (OutFilename_Info, argv[mark]);
                fp_out_info = fopen (OutFilename_Info, "wt");
                if (fp_out_info == NULL)
                {
                    fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", OutFilename_Info);
                    return FALSE;
                }
			}
			else if (fp_in_noisy)     /* Third argument string (input feature file) */
			{
				strcpy (InFilename_Fea, argv[mark]);
				fp_in_fea = fopen (InFilename_Fea, "rb");
				if (fp_in_fea == NULL)
				{
					fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", InFilename_Fea);
					return FALSE;
				}
			}
            else if (fp_in_clean)     /* Second argument string (input noisy wave file) */
            {
                strcpy (InFilename_Noisy, argv[mark]);
                fp_in_noisy = fopen (InFilename_Noisy, "rb");
                if (fp_in_noisy == NULL)
                {
                    fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", InFilename_Noisy);
                    return FALSE;
                }
            }
            else /* First argument string (input clean wave file) */
            {
                strcpy (InFilename, argv[mark]);
                fp_in_clean = fopen (argv[mark], "rb");
                if (fp_in_clean == NULL)
                {
                    fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", argv[mark]);
                    return FALSE;
                }
            }
            --argc;
            ++mark;
        }
    }

    if (!fp_in_clean || !fp_in_noisy || !fp_out_info || !fp_out_proc || !fp_in_fea )      /* Input and output files must be given */
    {
        fprintf (stderr, "ERROR:   Input and output files must be given!\r\n");
        return FALSE;
    }

    if (strcmp (InputKind, "NIST") &&
        strcmp (InputKind, "HTK") &&
        strcmp (InputKind, "RAW"))
    {
        fprintf (stderr, "ERROR:   Invalid input file format '%s'!\r\n", InputKind);
        return FALSE;
    }

    if (strcmp (InputKind, "RAW") && FsSpecified)
        fprintf (stderr, "WARNING:   Sampling frequency needs to be specified only for raw data files.\r\n");

    if (strcmp (InputKind, "RAW") && SwapSpecified)
        fprintf (stderr, "WARNING:   Byte swapping needs to be specified only for raw data files if necessary.\r\n");

	/*
    if ( Noc0 && NologE ) TextToParmKind( "MFCC", &OutParmKind );
    else if ( Noc0 ) TextToParmKind( "MFCC_E", &OutParmKind );
    else if ( NologE ) TextToParmKind( "MFCC_0", &OutParmKind );
	*/

    return TRUE;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: main
 *
 * PURPOSE:       Main front-end operations from input speech samples to output
 *                feature vectors. See embedded comments.
 *
 * INPUT:
 *   argc         Number of command line arguments (passed to ParseCommLine)
 *   argv         Array of command line arguments (passed to ParseCommLine)
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   TRUE         In case of any errors
 *   FALSE (0)    Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
extern int 
main (int argc, char *argv[])
{
    int i, j, k, TmpInt, CFBSize, CFBPointer;

    long FrameLength, FrameShift, FFTLength, FrameCounter = 0;

    float StartingFrequency, EnergyFloor_FB,
      EnergyFloor_logE, FloatBuffer[MAXWINDOWSIZE + 1], FloatBuffer_Noisy[MAXWINDOWSIZE+1],
      FloatWindow[MAXWINDOWSIZE / 2], *CircFloatBuffer, *CircFloatBuffer_Noisy;

	float temp_vectors[MAXWINDOWSIZE + 1];
	float smooth_buff[SMOOTH_WIN];
	float residual_noise_mean;
	float residual_noise_var;
	float residual_noise_max;
	float temp;

	float sum1, sum2, sum3;
	float clean_wave[MAXWINDOWSIZE + 1];

	float max_clean_spec;
	float max_noisy_spec;
	float max_denoise_spec;


    NIST_Header InNheader, InNheader_Noisy;
    HTK_Header InHheader, InHheader_Noisy;

    fprintf (stderr,"\r\nDSR Front-End v2.0\r\n");

    /*----------------*/
    /* Initialization */
    /*----------------*/
    EnergyFloor_FB = (float) exp ((double) ENERGYFLOOR_FB);
    EnergyFloor_logE = (float) exp ((double) ENERGYFLOOR_logE);
    InputByteOrder=NativeByteOrder;

    if (ParseCommLine (argc - 1, argv + 1))
    {
        /*-----------------------------------------------------------------*/
        /* Read input header, extract sampling frequency and byte ordering */
        /*-----------------------------------------------------------------*/
        if (!strcmp (InputKind, "NIST"))
        {
            if (!ReadNISTHeader (fp_in_clean, &InNheader))
            {
                fprintf (stderr, "ERROR:   Invalid NIST header !\r\n");
                goto _FaultExit;
            }
			if (!ReadNISTHeader (fp_in_noisy, &InNheader_Noisy))
			{
				fprintf (stderr, "ERROR:   Invalid NIST header !\r\n");
				goto _FaultExit;
			}
            SamplingFrequency = InNheader.SampleRate;
            if (strcmp (InNheader.SampleByteFormat, "10"))
                InputByteOrder = IEEE_LE;
            else
                InputByteOrder = IEEE_BE;
        }
        else if (!strcmp (InputKind, "HTK"))
        {
            if (!ReadHTKHeader (fp_in_clean, &InHheader, InputByteOrder!=IEEE_BE ))
            {
                fprintf (stderr, "ERROR:   Invalid HTK header !\r\n");
                goto _FaultExit;
            }
			if (!ReadHTKHeader (fp_in_noisy, &InHheader_Noisy, InputByteOrder!=IEEE_BE ))
			{
				fprintf (stderr, "ERROR:   Invalid HTK header !\r\n");
				goto _FaultExit;
			}
            /* 625->16kHz, 1250->8kHz, 909->11kHz */
            SamplingFrequency = (int) (10 * floor ((float) 1e6 /(float) InHheader.sampPeriod));
            InputByteOrder = IEEE_BE;
        }

		/* Read feature file */
		if (!ReadHTKHeader (fp_in_fea, &InHheader,0))
        //if (!ReadHTKHeader (fp_in_fea, &InHheader,1))
		{
            printf("读取HTK头信息出错!\n");
			fprintf (stderr, "ERROR:   Invalid HTK header !\r\n");
            goto _FaultExit;
        }

        /*-------------------------------------------------------------------*/
        /* Set parameters FrameLength, FrameShift and FFTLength according to */
        /* the current sampling frequency                                    */
        /*-------------------------------------------------------------------*/
        if (SamplingFrequency == SAMPLING_FREQ_1 * 1000)
        {
            FrameLength = FRAME_LENGTH_1;
            FrameShift = FRAME_SHIFT_1;
            FFTLength = FFT_LENGTH_1;
            StartingFrequency = STARTING_FREQ_1;
        }
        else if (SamplingFrequency == SAMPLING_FREQ_2 * 1000)
        {
            FrameLength = FRAME_LENGTH_2;
            FrameShift = FRAME_SHIFT_2;
            FFTLength = FFT_LENGTH_2;
            StartingFrequency = STARTING_FREQ_2;
        }
        else if (SamplingFrequency == SAMPLING_FREQ_3 * 1000)
        {
            FrameLength = FRAME_LENGTH_3;
            FrameShift = FRAME_SHIFT_3;
            FFTLength = FFT_LENGTH_3;
            StartingFrequency = STARTING_FREQ_3;
        }
        else
        {
            fprintf (stderr, "ERROR:   Invalid sampling frequency '%d'!\r\n",
		    SamplingFrequency);
            goto _FaultExit;
        }


		/* Load enhanced features */
		for(j=0; j<InHheader.nSamples; j++)
		{
			/* Load processed features */
			fread(temp_vectors, sizeof (float), (FFTLength/2+1),  fp_in_fea);
			for (i = 0; i <= FFTLength / 2; i++)
			{
				out_proc_vectors[i][j] = temp_vectors[i];
				if(out_proc_vectors[i][j] < ENERGYFLOOR_FB)
					out_proc_vectors[i][j] = EnergyFloor_FB;
				else
					out_proc_vectors[i][j] = (float) exp(out_proc_vectors[i][j]);
			}
		}

		/* Residual noise reduction*/
		if(SMOOTHPROCESS)
		{
			for (i = 0; i <= FFTLength / 2; i++)
			{
				residual_noise_mean = out_proc_vectors[i][0];
				residual_noise_var  = out_proc_vectors[i][0]*out_proc_vectors[i][0];
				residual_noise_max  = out_proc_vectors[i][0];
				for(j = 1; j < NOISE_FRAME_NUM ; j++)
				{
					residual_noise_mean += out_proc_vectors[i][j];
					residual_noise_var  += out_proc_vectors[i][j]*out_proc_vectors[i][j];
					if(residual_noise_max < out_proc_vectors[i][j])
						residual_noise_max = out_proc_vectors[i][j];
				}
				residual_noise_mean /= NOISE_FRAME_NUM;
				residual_noise_var = residual_noise_var/NOISE_FRAME_NUM - residual_noise_mean*residual_noise_mean;

				for(j = 0; j < SMOOTH_WIN; j++)
					smooth_buff[j] = out_proc_vectors[i][j];

				for(j = SMOOTH_WIN; j < (InHheader.nSamples - SMOOTH_WIN); j++)
				{
					if(out_proc_vectors[i][j]<residual_noise_max)
					{
						temp = out_proc_vectors[i][j];
						for(k = 1; k <= SMOOTH_WIN; k++)
						{
							if(temp>out_proc_vectors[i][j+k])
								temp = out_proc_vectors[i][j+k];

							if(temp>smooth_buff[k-1])
								temp = smooth_buff[k-1];

						}
						for(k = 1; k < SMOOTH_WIN; k++)
							smooth_buff[k-1] = smooth_buff[k];
						smooth_buff[SMOOTH_WIN-1] = out_proc_vectors[i][j];
						out_proc_vectors[i][j] = temp;
					}
					else
					{
						for(k = 1; k < SMOOTH_WIN; k++)
							smooth_buff[k-1] = smooth_buff[k];
						smooth_buff[SMOOTH_WIN-1] = out_proc_vectors[i][j];
					}
				}
			}

		}


	/*------------------------------------------------*/
	/* Memory allocation and initialization for input */
	/* circular float buffer                          */
	/*------------------------------------------------*/
	CFBSize=FrameLength;
	CircFloatBuffer=(float*)malloc(sizeof(float)*CFBSize);
	CircFloatBuffer_Noisy=(float*)malloc(sizeof(float)*CFBSize);
	if ( (!CircFloatBuffer) && (!CircFloatBuffer_Noisy) )
	{
	  fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
	  goto _FaultExit;
	}
	CFBPointer=0;

	segsnr = 0;
	segsnr_noisy = 0;
	lsd = 0;
	lsd_noisy = 0;

	/*-------------------------------------------------------*/
	/* Initialization of FE data structures and input buffer */
	/*-------------------------------------------------------*/
	InitializeHamming (FloatWindow, (int) FrameLength);

	ReadWave(fp_in_clean, CircFloatBuffer, CFBSize, CFBPointer, FrameLength-FrameShift,
		 (SwapSpecified || InputByteOrder != NativeByteOrder));

	ReadWave(fp_in_noisy, CircFloatBuffer_Noisy, CFBSize, CFBPointer, FrameLength-FrameShift,
		(SwapSpecified || InputByteOrder != NativeByteOrder));
	
	CFBPointer += (FrameLength-FrameShift);

	/*
	DCOffsetFilter( CircFloatBuffer, CFBSize, &CFBPointer, FrameLength-FrameShift );
	*/

        /*----------------------------------------------------------------*/
        /*                       Framing                                  */
        /*----------------------------------------------------------------*/
        while (ReadWave (fp_in_clean, CircFloatBuffer, CFBSize, CFBPointer%CFBSize, FrameShift,
			(SwapSpecified || InputByteOrder != NativeByteOrder))&&ReadWave (fp_in_noisy, CircFloatBuffer_Noisy, CFBSize, CFBPointer%CFBSize, FrameShift,
			(SwapSpecified || InputByteOrder != NativeByteOrder)))
        {
			CFBPointer = (CFBPointer+FrameShift)%CFBSize;
            FrameCounter++;

			sum1 = 0;
			sum2 = 0;
			for (i = 0; i < FrameLength; i++)
			{
                FloatBuffer[i] = CircFloatBuffer[(CFBPointer+i)%CFBSize];
				clean_wave[i]  = FloatBuffer[i];
				FloatBuffer_Noisy[i] = CircFloatBuffer_Noisy[(CFBPointer+i)%CFBSize];
				sum1 += FloatBuffer[i]*FloatBuffer[i];
				sum2 += (FloatBuffer_Noisy[i]-FloatBuffer[i])*(FloatBuffer_Noisy[i]-FloatBuffer[i]);
			}
			temp = 10*log10(sum1/sum2);

			if(temp>HIGHSEGSNR)
				temp=HIGHSEGSNR;
			if(temp<LOWSEGSNR)
				temp=LOWSEGSNR;
			segsnr_noisy += temp;

            /*-----------*/
            /* Windowing */
            /*-----------*/
            Window (FloatBuffer, FloatWindow, (int) FrameLength);
			Window (FloatBuffer_Noisy, FloatWindow, (int) FrameLength);

            /*-----*/
            /* FFT */
            /*-----*/

            /* Zero padding */
            for (i = FrameLength; i < FFTLength; i++)
			{
                FloatBuffer[i] = 0.0;
				FloatBuffer_Noisy[i] = 0.0;
			}

            /* Real valued, in-place split-radix FFT */
            TmpInt = (int) (log10 (FFTLength) / log10 (2));
            rfft (FloatBuffer, FFTLength, TmpInt); /*TmpInt = log2(FFTLength)*/
			rfft (FloatBuffer_Noisy, FFTLength, TmpInt); /*TmpInt = log2(FFTLength)*/

			/* Power spectrum */
			FloatBuffer[0] = (float) (FloatBuffer[0]*FloatBuffer[0]);  /* DC */
			for (i = 1; i < FFTLength / 2; i++)  /* pi/(N/2), 2pi/(N/2), ...,  (N/2-1)*pi/(N/2) */
				FloatBuffer[i] = (float) (FloatBuffer[i]*FloatBuffer[i]+FloatBuffer[FFTLength - i]*FloatBuffer[FFTLength - i]);
			FloatBuffer[FFTLength / 2] = (float) (FloatBuffer[FFTLength / 2]*FloatBuffer[FFTLength / 2]);  /* pi/2 */

			temp_vectors[0] = (float) (FloatBuffer_Noisy[0]*FloatBuffer_Noisy[0]);  /* DC */
			for (i = 1; i < FFTLength / 2; i++)  /* pi/(N/2), 2pi/(N/2), ...,  (N/2-1)*pi/(N/2) */
				temp_vectors[i] = (float) (FloatBuffer_Noisy[i]*FloatBuffer_Noisy[i]+FloatBuffer_Noisy[FFTLength - i]*FloatBuffer_Noisy[FFTLength - i]);
			temp_vectors[FFTLength / 2] = (float) (FloatBuffer_Noisy[FFTLength / 2]*FloatBuffer_Noisy[FFTLength / 2]);  /* pi/2 */


			for (i = 0; i <= FFTLength / 2; i++)  
			{
				clean_spec_vectors[i][FrameCounter-1] = FloatBuffer[i];
				noisy_spec_vectors[i][FrameCounter-1] = temp_vectors[i];
				denoise_spec_vectors[i][FrameCounter-1] = out_proc_vectors[i][FrameCounter-1];
			}

			/* Inverse process */
			/* post-processing */
			if(POSTPROCESS)
			{
				temp = FloatBuffer_Noisy[0]*FloatBuffer_Noisy[0];
				temp = log(temp);
				if(out_proc_vectors[0][FrameCounter-1]<(temp+THRESHOLD1))
					out_proc_vectors[0][FrameCounter-1] = temp+THRESHOLD1;
				for (i = 1; i <= SWITCHPOINT; i++) 
				{
					temp = FloatBuffer_Noisy[i]*FloatBuffer_Noisy[i]+FloatBuffer_Noisy[FFTLength - i]*FloatBuffer_Noisy[FFTLength - i];
					temp = log(temp);
					if(out_proc_vectors[i][FrameCounter-1]<(temp+THRESHOLD1))
						out_proc_vectors[i][FrameCounter-1] = temp+THRESHOLD1;
				}
				for (i = SWITCHPOINT+1; i < FFTLength/2; i++) 
				{
					temp = FloatBuffer_Noisy[i]*FloatBuffer_Noisy[i]+FloatBuffer_Noisy[FFTLength - i]*FloatBuffer_Noisy[FFTLength - i];
					temp = log(temp);
					if(out_proc_vectors[i][FrameCounter-1]<(temp+THRESHOLD2))
						out_proc_vectors[i][FrameCounter-1] = temp+THRESHOLD2;
				}
				temp = FloatBuffer_Noisy[FFTLength/2]*FloatBuffer_Noisy[FFTLength/2];
				temp = log(temp);
				if(out_proc_vectors[FFTLength/2][FrameCounter-1]<(temp+THRESHOLD2))
					out_proc_vectors[FFTLength/2][FrameCounter-1] = temp+THRESHOLD2;
			}
			

			/*Spectrum*/
			FloatBuffer_Noisy[0] = FloatBuffer_Noisy[0]*(sqrt(out_proc_vectors[0][FrameCounter-1])/fabs(FloatBuffer_Noisy[0]));
			FloatBuffer_Noisy[FFTLength/2] = FloatBuffer_Noisy[FFTLength/2]*(sqrt(out_proc_vectors[FFTLength/2][FrameCounter-1])/fabs(FloatBuffer_Noisy[FFTLength/2]));
			for (i = 1; i < FFTLength/2; i++) 
			{
				temp = (float) sqrt((FloatBuffer_Noisy[i]*FloatBuffer_Noisy[i]+FloatBuffer_Noisy[FFTLength - i]*FloatBuffer_Noisy[FFTLength - i]));
				temp = (float) (sqrt(out_proc_vectors[i][FrameCounter-1])/temp);
				FloatBuffer_Noisy[i] *= temp;
				FloatBuffer_Noisy[FFTLength - i] *= temp;
			}

			rifft(FloatBuffer_Noisy, FFTLength, TmpInt); 

			for (i = 0; i < FrameLength; i++)  
			{
				temp_vectors[i] = FloatBuffer_Noisy[i];
			}
			DeWindow (temp_vectors, FloatWindow, (int) FrameLength);
			sum2 = 0;
			for (i = 0; i < FrameLength; i++)
			{
				sum2 += (temp_vectors[i]-clean_wave[i])*(temp_vectors[i]-clean_wave[i]);
			}
			temp = 10*log10(sum1/sum2);

			if(temp>HIGHSEGSNR)
				temp=HIGHSEGSNR;
			if(temp<LOWSEGSNR)
				temp=LOWSEGSNR;
			segsnr += temp;

			if(OLA_KIND==1)
				Window (FloatBuffer_Noisy, FloatWindow, (int) FrameLength);
			else
				DeWindow (FloatBuffer_Noisy, FloatWindow, (int) FrameLength);

			/*Store*/
			for (i = 0; i < FrameLength; i++)
			{
				out_proc_vectors[i][FrameCounter-1] = FloatBuffer_Noisy[i];
			}
            
            /*---------------------------*/
            /* Display processing status */
            /*---------------------------*/
            if (!QuietMode && !(FrameCounter % PRINTMOD))
            {
                fprintf (stderr, "\rProcessing status: %ld frames ...", FrameCounter);
                fflush (stderr);
            }

        }

	/*LSD*/
	max_clean_spec = clean_spec_vectors[0][0];
	max_noisy_spec = noisy_spec_vectors[0][0];
	max_denoise_spec = denoise_spec_vectors[0][0];

	for(i=0; i< FrameCounter; i++)
	{
		for (j = 0; j <= FFTLength / 2; j++)  
		{
			if(max_clean_spec<clean_spec_vectors[j][i])
				max_clean_spec = clean_spec_vectors[j][i];

			if(max_noisy_spec<noisy_spec_vectors[j][i])
				max_noisy_spec = noisy_spec_vectors[j][i];

			if(max_denoise_spec<denoise_spec_vectors[j][i])
				max_denoise_spec = denoise_spec_vectors[j][i];
		}
	}

	max_clean_spec *= pow(10.0, -DYNRANGE/10.0);
	max_noisy_spec *= pow(10.0, -DYNRANGE/10.0);
	max_denoise_spec *= pow(10.0, -DYNRANGE/10.0);


	for(i=0; i< FrameCounter; i++)
	{
		for (j = 0; j <= FFTLength / 2; j++)  
		{
			if(max_clean_spec>clean_spec_vectors[j][i]) 
				clean_spec_vectors[j][i] = max_clean_spec;

			if(max_noisy_spec>noisy_spec_vectors[j][i])
				noisy_spec_vectors[j][i] = max_noisy_spec;

			if(max_denoise_spec>denoise_spec_vectors[j][i])
				denoise_spec_vectors[j][i] = max_denoise_spec;
		}
	}

	for(i=0; i< FrameCounter; i++)
	{
		sum3 = 0;
		for (j = 0; j <= FFTLength / 2; j++)  
		{
			temp = 10*log10(noisy_spec_vectors[j][i]/clean_spec_vectors[j][i]);
			sum3 += temp*temp;
		}
		sum3 /= (FFTLength/2+1);
		lsd_noisy += sqrt(sum3);
	}

	for(i=0; i< FrameCounter; i++)
	{
		sum3 = 0;
		for (j = 0; j <= FFTLength / 2; j++)  
		{
			temp = 10*log10(denoise_spec_vectors[j][i]/clean_spec_vectors[j][i]);
			sum3 += temp*temp;
		}
		sum3 /= (FFTLength/2+1);
		lsd += sqrt(sum3);
	}


	/*Waveform Reconstruction*/
	waveform_num = FrameCounter*FrameShift+(FrameLength-FrameShift);
	for(i=0; i<waveform_num; i++)
	{
		waveform_count[i] = 0.0;
		out_proc_wave[i]  = 0.0;
	}
	
	for(i=0; i<FrameCounter; i++)
	{
		for(j=0; j<FrameLength; j++)
		{
			if(OLA_KIND==1)
			{
				if(j<FrameLength/2)
					waveform_count[i*FrameShift+j] += FloatWindow[j]*FloatWindow[j];
				else
					waveform_count[i*FrameShift+j] += FloatWindow[FrameLength - 1 -j]*FloatWindow[FrameLength - 1 -j];
			}
			else
			{
				waveform_count[i*FrameShift+j] += 1.0;
			}
			out_proc_wave[i*FrameShift+j] += out_proc_vectors[j][i];
		}
	}

	for(i=0; i<waveform_num; i++)
	{
		out_proc_wave[i] /= waveform_count[i];
		waveform[i] = (short) out_proc_wave[i];
	}
	fwrite(waveform, sizeof(short), waveform_num, fp_out_proc);

	/* Output Information */
	segsnr /= FrameCounter;
	lsd /= FrameCounter;
	segsnr_noisy /= FrameCounter;
	lsd_noisy /= FrameCounter;

	fprintf(fp_out_info,"Segmental SNR:\n");
	fprintf(fp_out_info,"%f\n",segsnr);
	fprintf(fp_out_info,"Log-Spectral Distortion:\n");
	fprintf(fp_out_info,"%f\n",lsd);

	if(NoisyInfo)
	{
		strcpy(OutFilename_Noisy_Info,InFilename_Noisy);
		strcat(OutFilename_Noisy_Info,".info");

		fp_out_noisy_info = fopen (OutFilename_Noisy_Info, "wt");
		if (fp_out_noisy_info == NULL)
		{
			fprintf (stderr, "ERROR:   Could not open file '%s' !\r\n", OutFilename_Noisy_Info);
			return FALSE;
		}

		fprintf(fp_out_noisy_info,"Segmental SNR:\n");
		fprintf(fp_out_noisy_info,"%f\n",segsnr_noisy);
		fprintf(fp_out_noisy_info,"Log-Spectral Distortion:\n");
		fprintf(fp_out_noisy_info,"%f\n",lsd_noisy);

		fclose(fp_out_noisy_info);
	}
	

	/*----------------*/
	/* Memory release */
	/*----------------*/
	free(CircFloatBuffer);
	free(CircFloatBuffer_Noisy);

        /*------------------------------*/
        /* Close input and output files */
        /*------------------------------*/
        fclose (fp_in_clean);
		fclose (fp_in_noisy);
		fclose (fp_in_fea);
        fclose (fp_out_info);
		fclose (fp_out_proc);

        /*----------------------*/
        /* Display final status */
        /*----------------------*/
        if (!QuietMode)
            fprintf (stderr, "\rProcessed: %ld Frames.                      \r\n", FrameCounter);

    }                           /* if ( ParsCommLine... ) */
    else
    {
        fprintf (stderr, "\r\nUSAGE:");
        fprintf (stderr, "   %s infile infile_feature outfile_ori_wav outfile_proc_wav [options]\r\n", argv[0]);
        fprintf (stderr, "\r\nOPTIONS:\r\n");
        fprintf (stderr, "     -q            Quiet Mode                                 (%s)\r\n", QuietMode ? "TRUE" : "FALSE");
        fprintf (stderr, "     -F    format  Input file format (NIST,HTK,RAW)           (%s)\r\n", InputKind);
        fprintf (stderr, "     -fs   freq    Sampling frequency in kHz (%d,%d,%d)        (%d)\r\n", SAMPLING_FREQ_1, SAMPLING_FREQ_2, SAMPLING_FREQ_3, SamplingFrequency / 1000);
        fprintf (stderr, "     -swap         Change input byte ordering                 (%s)\r\n", SwapSpecified ? "Swapped" : "Native");
        fprintf (stderr, "                   (Native byte ordering is %s)\r\n", NativeByteOrder ? "ieee-be" : "ieee-le");
		fprintf (stderr, "     -noh          No HTK header to output file               (%s)\r\n", NoOutHeaderSpecified ? "TRUE" : "FALSE" );
		fprintf (stderr, "     -noc0         No c0 coefficient to output feature vector (%s)\r\n", Noc0 ? "TRUE" : "FALSE" );
		fprintf (stderr, "     -nologE       No logE component to output feature vector (%s)\r\n", NologE ? "TRUE" : "FALSE" );

    }

    fprintf (stderr, " Copyright (c) 1998 Nokia Research Center, Tampere, Finland. All rights reserved.\r\n");
    return FALSE;

  _FaultExit:
    if (fp_in_clean)
        fclose (fp_in_clean);
	if (fp_in_noisy)
		fclose (fp_in_noisy);
	if (fp_in_fea)
        fclose (fp_in_fea);
    if (fp_out_info)
        if (fclose (fp_out_info))
            fprintf (stderr, "Can't rewrite output file!\r\n");
	if (fp_out_proc)
        if (fclose (fp_out_proc))
            fprintf (stderr, "Can't rewrite output file!\r\n");
    return TRUE;
}





