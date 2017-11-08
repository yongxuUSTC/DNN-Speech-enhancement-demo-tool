/*----------------------------------------------------------------------------
 *
 * Mel Cepstrum Parametrisation (Track 1 Front End) v2.0 for Distributed
 * Speech Recognition (DSR), file I/O functions
 *
 * FILE NAME: fileio.c
 * PURPOSE:   General purpose function package for DSR file I/O operations.
 *            Various input file formats (RAW, NIST, HTK) are supported.
 *            Output file format is always HTK. Package consists of functions
 *            for reading different file headers and input data and writing
 *            output HTK header and output data.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"

/*---------------------------
 * Local Function Prototypes
 *---------------------------*/
int ParseHTKParmKind (int parmKind, char *outstr);
void Swap16 ( short *Short );
void Swap32 ( long *Long );

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ReadNISTHeader
 *
 * PURPOSE:       Reads NIST header of an input file and stores header
 *                information in output data structure. Called by main
 *                function.
 *
 * INPUT:
 *   fp_in        Input waveform file handle
 *   header       Pointer to memory segment allocated for output data structure
 *
 * OUTPUT
 *                Header information stored in output data structure pointed
 *                to by *header*
 *
 * RETURN VALUE
 *   0            In case of any errors
 *   1            Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
int 
ReadNISTHeader (FILE * fp_in, NIST_Header * header)
{
    char str[199], s1[30], s2[30], s3[30];

    fseek (fp_in, 0L, SEEK_SET);
    if (!fgets (str, sizeof (str), fp_in))
        return 0;               /* Empty file */

    str[sizeof (header->Type) - 1] = '\0';
    strcpy (header->Type, str);
    str[4] = '\0';
    if (strcmp (str, "NIST"))
        return 0;               /* No NIST stamp */

    if (!fgets (str, sizeof (str), fp_in))
        return 0;               /* Unexpected end of file */

    header->HeaderSize = atoi (str);

    while (fgets (str, sizeof (str), fp_in))
    {
        sscanf (str, "%s %s %s", s1, s2, s3);
        if (!strcmp (s1, "end_head"))
        {
            fseek (fp_in, (long) (header->HeaderSize), SEEK_SET);
            return 1;
        }
        else if (!strcmp (s1, "database_id"))
            strcpy (header->DatabaseID, s3);
        else if (!strcmp (s1, "database_version"))
            strcpy (header->DatabaseVersion, s3);
        else if (!strcmp (s1, "utterance_id"))
            strcpy (header->UtteranceID, s3);
        else if (!strcmp (s1, "channel_count"))
            header->ChannelCount = atoi (s3);
        else if (!strcmp (s1, "sample_count"))
            header->SampleCount = atol (s3);
        else if (!strcmp (s1, "sample_rate"))
            header->SampleRate = atol (s3);
        else if (!strcmp (s1, "sample_min"))
            header->SampleMin = atoi (s3);
        else if (!strcmp (s1, "sample_max"))
            header->SampleMax = atoi (s3);
        else if (!strcmp (s1, "sample_n_bytes"))
            header->SamplenBytes = atoi (s3);
        else if (!strcmp (s1, "sample_byte_format"))
            strcpy (header->SampleByteFormat, s3);
        else if (!strcmp (s1, "sample_sig_bits"))
            header->SampleSigBits = atoi (s3);
        else if (!strcmp (s1, "sample_coding"))
            strcpy (header->SampleCoding, s3);
        else if (!strcmp (s1, "sample_checksum"))
            header->SampleChecksum = atol (s3);
    }
    return 0;                   /* Incorrect header format */
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ReadHTKHeader
 *
 * PURPOSE:       Reads HTK header of an input file and stores header
 *                information in output data structure. Called by main
 *                function.
 *
 * INPUT:
 *   fp_in        Input waveform file handle
 *   header       Pointer to memory segment allocated for output data structure
 *   swap         Nonzero if byte swapping is needed
 *
 * OUTPUT
 *                Header information stored in output data structure pointed
 *                to by *header*
 *
 * RETURN VALUE
 *   0            In case of any errors
 *   1            Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
int 
ReadHTKHeader (FILE * fp_in, HTK_Header * header, int Swap )
{

    char str[20];

    fseek (fp_in, 0L, SEEK_SET);
    if (!fread (&(header->nSamples), sizeof (long), 1, fp_in))
	{  printf("第一个错误\n") ;return 0;}
    if (!fread (&(header->sampPeriod), sizeof (long), 1, fp_in))
	{  printf("第二个错误\n") ;return 0;}
    if (!fread (&(header->sampSize), sizeof (short), 1, fp_in))
          {  printf("第三个错误\n") ;return 0;}
    if (!fread (&(header->sampKind), sizeof (short), 1, fp_in))
         {  printf("第四个错误\n") ;return 0;}

    if ( Swap )
    {
        Swap32 (&(header->nSamples));
        Swap32 (&(header->sampPeriod));
        Swap16 (&(header->sampSize));
        Swap16 (&(header->sampKind));
    }
    
//    if (header->nSamples < 0 || header->sampPeriod < 0 ||
//        header->sampSize < 0 ||
//        !ParseHTKParmKind (header->sampKind, str))
//       {  printf("第五个错误\n") ;return 0;}

	//printf("header->nSamples=%d\n",header->nSamples) ;
	//printf("header->sampSize=%d\n",header->sampSize) ;
	//printf("header->sampKind=%d\n",header->sampKind) ;
	//printf("header->sampPeriod=%d\n",header->sampPeriod) ;

	if(header->nSamples < 0)
	{  printf("header->nSamples=%d,第五个错误\n",header->nSamples) ;return 0;}
	
	else if(header->sampSize < 0)
		{  printf("header->sampSize=%d,第7个错误\n",header->sampSize) ;return 0;}
	
		else if(!ParseHTKParmKind (header->sampKind, str))
		{  printf("header->sampKind=%d,第8个错误\n",header->sampKind) ;return 0;}

	else if (header->sampPeriod < 0)
		{  printf("header->sampPeriod=%d,第6个错误\n",header->sampPeriod) ;return 0;}


    else
        return 1;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: WriteHTKHeader
 *
 * PURPOSE:       Writes HTK header stored in input data structure into an
 *                output file. Called by main function.
 *
 * INPUT:
 *   fp_out       Output feature file handle
 *   header       Pointer to HTK data structure
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
WriteHTKHeader (FILE * fp_out, HTK_Header * header)
{
    fseek (fp_out, 0L, SEEK_SET);
    fwrite (&(header->nSamples), sizeof (long), 1, fp_out);
    fwrite (&(header->sampPeriod), sizeof (long), 1, fp_out);
    fwrite (&(header->sampSize), sizeof (short), 1, fp_out);
    fwrite (&(header->sampKind), sizeof (short), 1, fp_out);
    fseek (fp_out, 0L, SEEK_END);
}

/*---------------------------------------------------------------------------
 * FUNCTION NAME: WriteHTKFeature
 *
 * PURPOSE:       Writes one feature vector into output file.
 *                Called by main function.
 *
 * INPUT:
 *   fp_out       Output feature file handle
 *   out          Pointer to feature vector
 *   fea_len      Feature vector length
 *
 * OUTPUT
 *   none
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void 
WriteHTKFeature (FILE * fp_out, float *out, short fea_len)
{
    fwrite (out, sizeof (float), fea_len, fp_out);
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ReadWave
 *
 * PURPOSE:       Reads waveform into memory
 *
 * INPUT:
 *   fp_in        Input waveform file handle
 *   CircBuff     Circular float buffer
 *   BSize        Buffer size
 *   BPointer     Buffer pointer
 *   nSamples     Number of samples to be read
 *   Swap         Nonzero if byte swapping is necessary
 *
 * OUTPUT
 *                Waveform stored in circular float buffer. Buffer pointer
 *                remains untouched (points to the first sample read)
 *
 * RETURN VALUE
 *   0            In case of any errors
 *   1            Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
int
ReadWave (FILE *fp_in, float *CircBuff, int BSize, int BPointer, int nSamples, int Swap )
{
  short s;
  int i;

  for ( i=0; i<nSamples; i++ )
    {
    if (fread (&s, sizeof (short), 1, fp_in) != 1) return 0;
    if (Swap) s = ((s & 0x00ff) << 8) | ((s & 0xff00) >> 8);
    CircBuff[BPointer]=(float)s;
    BPointer=(BPointer+1)%BSize;
    }
  return 1;
}



/*-----------------
 * Local Functions
 *-----------------*/

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ParseHTKParmKind
 *
 * PURPOSE:       Parses HTK parameter kind integer and converts in to ASCII
 *                format.
 *
 * INPUT:
 *   parmKind     Input parameter kind integer from HTK header
 *   outstr       Pointer to ASCII buffer
 *
 * OUTPUT
 *                ASCII parameter kind stored in *outstr*
 *
 * RETURN VALUE
 *   0            If input parameter kind in incorrect
 *   1            Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
int 
ParseHTKParmKind (int parmKind, char *outstr)
{
    int tmp;

    tmp = parmKind & 0x003f;    /* Basic Parameter Kind Codes */
    switch (tmp)
    {
    case 0:
        strcpy (outstr, "WAVEFORM");
        break;
    case 1:
        strcpy (outstr, "LPC");
        break;
    case 2:
        strcpy (outstr, "LPREFC");
        break;
    case 3:
        strcpy (outstr, "LPCEPSTRA");
        break;
    case 4:
        strcpy (outstr, "LPDELCEP");
        break;
    case 5:
        strcpy (outstr, "IREFC");
        break;
    case 6:
        strcpy (outstr, "MFCC");
        break;
    case 7:
        strcpy (outstr, "FBANK");
        break;
    case 8:
        strcpy (outstr, "MELSPEC");
        break;
    case 9:
        strcpy (outstr, "USER");
        break;
    case 10:
        strcpy (outstr, "DISCRETE");
        break;
    default:
        return 0;
    }
    tmp = parmKind & 0xffc0;    /* Bit-encoded Qualifiers */
    if (tmp & 000100)
        strcat (outstr, "_E");
    if (tmp & 000200)
        strcat (outstr, "_N");
    if (tmp & 000400)
        strcat (outstr, "_D");
    if (tmp & 001000)
        strcat (outstr, "_A");
    if (tmp & 002000)
        strcat (outstr, "_C");
    if (tmp & 004000)
        strcat (outstr, "_Z");
    if (tmp & 010000)
        strcat (outstr, "_K");
    if (tmp & 020000)
        strcat (outstr, "_0");
    return 1;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: TextToParmKind
 *
 * PURPOSE:       Converts ASCII HTK parameter description to parameter kind
 *                integer to be written into HTK header
 *
 * INPUT:
 *   instr        Pointer to input parameter string
 *   parmKind     Pointer to output parameter kind integer
 *
 * OUTPUT
 *                Output parameter kind integer
 *
 * RETURN VALUE
 *   0            If input parameter string is incorrect
 *   1            Otherwise
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
int
TextToParmKind( char *instr, int *parmKind )
{
 char   workstr[40];
 int    idx;

 *parmKind=0;
 strcpy( workstr, instr );
 idx= (int) strlen(workstr);

 while ( workstr[idx-2]=='_' )
   {
   workstr[idx]='\0';
   switch ( workstr[idx-1] )
     {
     case 'E' : *parmKind+=000100;
                break;
     case 'N' : *parmKind+=000200;
                break;
     case 'D' : *parmKind+=000400;
                break;
     case 'A' : *parmKind+=001000;
                break;
     case 'C' : *parmKind+=002000;
                break;
     case 'Z' : *parmKind+=004000;
                break;
     case 'K' : *parmKind+=010000;
                break;
     case '0' : *parmKind+=020000;
                break;
     default  : return 0;
     } 
   idx-=2;
   if ( idx<2 ) return 0;
   }
   workstr[idx]='\0';

   if ( !strcmp( workstr, "WAVEFORM" ) ) ;
   else if ( !strcmp( workstr, "LPC" ) ) *parmKind+=1;
   else if ( !strcmp( workstr, "LPREFC" ) ) *parmKind+=2;
   else if ( !strcmp( workstr, "LPCEPSTRA" ) ) *parmKind+=3;
   else if ( !strcmp( workstr, "LPDELCEP" ) ) *parmKind+=4;
   else if ( !strcmp( workstr, "IREFC" ) ) *parmKind+=5;
   else if ( !strcmp( workstr, "MFCC" ) ) *parmKind+=6;
   else if ( !strcmp( workstr, "FBANK" ) ) *parmKind+=7;
   else if ( !strcmp( workstr, "MELSPEC" ) ) *parmKind+=8;
   else if ( !strcmp( workstr, "USER" ) ) *parmKind+=9;
   else if ( !strcmp( workstr, "DISCRETE" ) ) *parmKind+=10;
   else return 0;

   return 1;
}


/*----------------------------------------------------------------------------
 * FUNCTION NAME: Swap16
 *
 * PURPOSE:       Swaps two bytes
 *
 * INPUT:
 *   Short        Pointer to input and output short variable
 *
 * OUTPUT
 *                Byte swapped short variable pointed by *Short*
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void Swap16 ( short *Short )
{
  *Short = ((*Short&0x0000ff00L)<<8 )|((*Short&0x00ff0000L)>>8 );
}


/*----------------------------------------------------------------------------
 * FUNCTION NAME: Swap32
 *
 * PURPOSE:       Reorders four bytes
 *
 * INPUT:
 *   Long         Pointer to input and output long variable
 *
 * OUTPUT
 *                Byte swapped long variable pointed by *Long*
 *
 * RETURN VALUE
 *   none
 *
 * Copyright (c) 1998 Nokia Research Center, Tampere, Finland
 *---------------------------------------------------------------------------*/
void Swap32 ( long* Long )
{
  *Long = ((*Long&0x000000ffL)<<24 )| \
          ((*Long&0x0000ff00L)<<8 )| \
          ((*Long&0x00ff0000L)>>8 )| \
          ((*Long&0xff000000L)>>24 );
}
