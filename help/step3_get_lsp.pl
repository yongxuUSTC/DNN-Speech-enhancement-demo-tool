use warnings;

$scp="tmp\\yqzu_raw_16k_noisy.scp";
open(IN1,$scp);
while(<IN1>)
{
	chomp;
	$line=$_;
	if($line=~m/\.raw/)
	{
		$inline=$line;
		$outline=$line;
		$outline=~s/\.raw/\.lsp/;
	system("help\\Wav2LogSpec.exe -F RAW -fs 16 $inline $outline");
		}

	#die;
	}
