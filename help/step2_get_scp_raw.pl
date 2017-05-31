use warnings;

#my ($path, $oscp) = @ARGV;
$path="wav_lsp";
$oscp="tmp\\yqzu_raw_16k_noisy.scp";

open(OUT1,">$oscp");
opendir(DIR,$path);
foreach $f (readdir DIR)
{
	if($f=~m/\.raw/)
	{
		print OUT1 "$path\\$f\n";
		}
	}
