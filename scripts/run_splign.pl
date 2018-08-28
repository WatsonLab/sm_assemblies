#!/usr/bin/perl

use Bio::SearchIO;

my %int;
# interrupted file
open(INT, $ARGV[0]);
while(<INT>) {
	chomp();
	$int{$_}++;
}
close INT;

my $dir = $ARGV[3];

# top level dirs
mkdir($dir) unless (-d $dir);
mkdir("$dir/CDS") unless (-d "$dir/CDS");

# BLAST file
my $in = Bio::SearchIO->new(-file => $ARGV[1], -format => 'blast');
while(my $result = $in->next_result) {
	
	# get top two hits
	my $hit = $result->next_hit;
	next unless defined $hit;
	my $hit2 = $result->next_hit;

	# get query details	
	my $qv = $result->query_name;
	my $qn = $qv;
	$qn =~ s/\.\d+//;

	next unless (exists $int{$qn});

	my @hits = ();
        push(@hits, $hit) if defined ($hit);
        push(@hits, $hit2) if defined ($hit2);

	foreach $hit (@hits) {

                my $hn = $hit->name;

                # CDS
                unless (-d "$dir/CDS/$qv") {
                        mkdir("$dir/CDS/$qv");
                }

                system("samtools faidx Homo_sapiens.GRCh38.cds.all.fa.gz \"$qv\" > $dir/CDS/$qv/$qv.fa");
                # the Assembly
                system("samtools faidx $ARGV[2] \"$hn\" > $dir/CDS/$qv/$hn.fa");
                print STDERR "Trying splign on $qv and CDS $hn\n";
                system("./splign -type mrna -query $dir/CDS/$qv/$qv.fa -subj $dir/CDS/$qv/$hn.fa -aln $dir/CDS/$qv/$qv.$hn.aln -log $dir/CDS/$qv/$qv.$hn.splign.log > $dir/CDS/$qv/$qv.hn.out");
                system("rm $dir/CDS/$qv/$hn.fa");

                print "$qv\t$qn\t$hn\n";
        }
}

