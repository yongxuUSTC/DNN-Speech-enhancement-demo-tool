use warnings;

#my ($path, $oscp) = @ARGV;
$path="wav_lsp";
$oscp="tmp\\yqzu_lsp_16k_noisy.scp";

open(OUT1,">$oscp");
opendir(DIR,$path);
foreach $f (readdir DIR)
{
	if($f=~m/\.lsp/)
	{
		print OUT1 "$path\\$f\n";
		}
	}
system("copy tmp\\yqzu_lsp_16k_noisy.scp scp\\yqzu_lsp_16k_noisy.scp");
