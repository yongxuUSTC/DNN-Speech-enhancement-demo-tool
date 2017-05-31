use warnings;

##noisy
$inpath="wav_lsp";
$outpath="wav_lsp";

if (!(-e $outpath)) {system("mkdir $outpath");}

#$inpath="D\:\\user\\jundu\\xuyong\\yanping_10h\\raw_8k";
#$outpath="D\:\\user\\jundu\\xuyong\\yanping_10h\\lsp_8k";


#open(OUT1,">scp\\yanping_all_lsp_babble_8k_N2SNR10.scp");
#open(OUT1,">scp\\yanping_all_lsp_8k_clean.scp");
opendir(DIR,$inpath);
foreach $f (readdir DIR)
{
	if($f=~m/^\./){next;}
	$outf=$f;
	$outf=~s/\.wav/\.raw/;
	
	system("help\\WAV2RAW.exe $inpath\\$f $outpath\\$outf");
#	print OUT1 "$outpath\\$outf\n";
}

###clean
#$inpath="D:\\user\\jundu\\xuyong\\cmu_one_speaker\\CMU_ARCTIC\\cmu_us_slt_arctic-0.95-release\\cmu_us_slt_arctic\\raw_8k";
#$outpath="D:\\user\\jundu\\xuyong\\cmu_one_speaker\\CMU_ARCTIC\\cmu_us_slt_arctic-0.95-release\\cmu_us_slt_arctic\\lsp_8k";
#
#open(OUT1,">cmu_arctic_female_lsp_clean_8k.scp");
#opendir(DIR,$inpath);
#foreach $f (readdir DIR)
#{
#	if($f=~m/^\./){next;}
#	$outf=$f;
#	$outf=~s/\.raw/\.lsp/;
#	
#	system("D\:\\user\\jundu\\AE_SE\\Bin\\Wav2LogSpec.exe -F RAW -fs 8 $inpath\\$f $outpath\\$outf");
#	print OUT1 "$outpath\\$outf\n";
#}