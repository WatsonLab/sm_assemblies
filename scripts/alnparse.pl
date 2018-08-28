#!/usr/bin/perl

use Data::Dumper;
use Bio::Perl;

my $fullout = undef;
my $fullout = $ARGV[1];

my @hits;
my $h = undef;
my $e = 0;
open(IN, $ARGV[0]);

my @entirefile = <IN>;

for ($i=0;$i<@entirefile;$i++) {

	my $line = $entirefile[$i];

	chomp($line);

	if ($line =~ m/^>/) {
		# we have a new hit!
		#print "there is a hit!\n$line\n";		
		# if we have a previous hit, store it
		# and undef the variable $h
		push(@hits, &augment($h)) if defined ($h);			
		$h = undef;	
 	
		# reset exon counter to 0
		$e = 0;

		if ($line =~ m/^>(\S+)\s+(\S+)\s+(\S+)/) {
			$h->{details}->{hitstrand} = $1;
			$h->{details}->{queryname} = $2;
			$h->{details}->{hitname}   = $3;
		}
	}

	if ($line =~ m/ Exon (\d+) \((\d+)-(\d+),(\d+)-(\d+)\) Len = (\d+) Identity = (\d+)/) {
	
		$e = $1;

		#print "Found exon $e\n";
		$h->{exons}->{$e}->{qstart} = $2;
		$h->{exons}->{$e}->{qend}   = $3;
		$h->{exons}->{$e}->{hstart} = $4;
		$h->{exons}->{$e}->{hend}   = $5;
		$h->{exons}->{$e}->{alen}   = $6;
		$h->{exons}->{$e}->{fid}    = $7;

		while ($entirefile[$i+1] !~ m/ Exon/ && $entirefile[$i+1] !~ m/^>/ && $entirefile[$i+2] !~ m/^>/) {
			#print "while $i\n";
			my $prot = $entirefile[++$i];
			my $dna1 = $entirefile[++$i];
			my $mat  = $entirefile[++$i];
			my $dna2 = $entirefile[++$i];
			my $blan = $entirefile[++$i];
			
			last unless ($dna1 =~ m/[AGCT]/);
	
			chomp($prot);
			chomp($dna1);
			chomp($mat);
			chomp($dna2);
					
			# trim off left-most detail
			my $ltrim_string;
			if ($dna1 =~ m/^(\s+\d+\s+)/) {
				$ltrim_string = $1;
			}

			my $ltrim = length $ltrim_string;

			$prot = substr($prot,$ltrim);
			$dna1 = substr($dna1,$ltrim);
			$mat  = substr($mat, $ltrim);
			$dna2 = substr($dna2,$ltrim);

			# trim off left-most dots
			my $dltrim_string;
                        if ($dna1 =~ m/^(\.+)/) {
                                $dltrim_string = $1;
                        }

                        my $dltrim = length $dltrim_string;

                        $prot = substr($prot,$dltrim);
                        $dna1 = substr($dna1,$dltrim);
                        $mat  = substr($mat, $dltrim);
                        $dna2 = substr($dna2,$dltrim);

			# trim off right most dots
			my $drtrim_string;
                        if ($dna1 =~ m/(\.+)$/) {
                                $drtrim_string = $1;
                        }

                        my $drtrim = length $drtrim_string;

                        $prot = substr($prot, 0, -1 * $drtrim) if ($drtrim > 0);
                        $dna1 = substr($dna1, 0, -1 * $drtrim) if ($drtrim > 0);
                        $mat  = substr($mat,  0, -1 * $drtrim) if ($drtrim > 0);
                        $dna2 = substr($dna2, 0, -1 * $drtrim) if ($drtrim > 0);

			if ($prot =~ m/[A-Z]/) {
				$h->{hasprot} = 1;
			} else {
				$h->{hasprot} = 0;
			}	

			$h->{prot} .= $prot;
			$h->{dna1} .= $dna1;
			$h->{mat}  .= $mat;
			$h->{dna2} .= $dna2;

			last unless ($dna1 =~ m/[AGCT]/);

		}	
	}
}

push(@hits, &augment($h)) if defined ($h);


# if there is only one hit, process it
if (scalar(@hits) == 1) {
	&process_hits(@hits);
	exit(0);
}

# when we get to here, there are multiple hits

# count proteins
my @prothits = ();
foreach $h (@hits) {
	if ($h->{hasprot} == 1) {
		push(@prothits, $h);
	}
}

exit unless (@prothits > 0);

# if there is only one protein hit
# process it
if (scalar(@prothits)==1) {
	&process_hits(@prothits);
	exit(0);
}

# when we get to here, there are
# multiple hits with proteins

# sort on # of indels
my @indel_sorted = sort {($a->{dels}+$a->{ins}) <=> ($b->{dels}+$b->{ins})} @prothits;
my @min_indels;
push(@min_indels, $indel_sorted[0]);
for ($i=1;$i<@indel_sorted;$i++) {
	if (($indel_sorted[$i]->{dels}+$indel_sorted[$i]->{ins}) <= ($indel_sorted[0]->{dels}+$indel_sorted[0]->{ins})) {
		push(@min_indels, $indel_sorted[$i]);
	}
}

# process if there is a clear winner
if (scalar(@min_indels)==1) {
	&process_hits(@min_indels);
	exit(0);
}

# if we get here, there are multiple
# protein hits and there are at
# least two with the same minimum # of indels

# sort on # of mismatches
my @match_sorted = sort {$a->{mis} <=> $b->{mis}} @prothits;

my @min_match;
push(@min_match, $match_sorted[0]);
for ($i=1;$i<@match_sorted;$i++) {
        if ($match_sorted[$i]->{mis} <= $match_sorted[0]->{mis}) {
                push(@min_match, $match_sorted[$i]);
        }
}

# process if there is a clear winner
if (scalar(@min_match)==1) {
        &process_hits(@min_match);
        exit(0);
} 


sub process_hits {

	my @hits = @_;

	foreach $h (@hits) {

        	next unless ( $h->{hasprot} == 1);

		print $h->{details}->{queryname}, "\t", $h->{details}->{hitname}, "\t";
		print $h->{exons}->{1}->{qstart}, "\t", length($h->{cleandna1}), "\t";
		print $h->{mis}, "\t", $h->{misb}, "\t";
        	print $h->{ins}, "\t", $h->{insb}, "\t";
        	print $h->{dels}, "\t", $h->{delb}, "\t";
        
		my $pstring = translate_as_string($h->{cleandna2});
		$pstring =~ s/\*/__*__/g;
		#print $pstring, "\n";
		
		print $pstring, "\n";

		if (defined $fullout) {
                        print $h->{prot}, "\n";
                        print $h->{dna1}, "\n";
                        print $h->{mat},  "\n";
                        print $h->{dna2}, "\n";
                }

	}
}

sub augment {

	my $h = shift;

	my $dna1 = $h->{dna1};
        my $mat  = $h->{mat};
        my $dna2 = $h->{dna2};

	# dashed in dna1 indicate insertions in reference
	my $ins  = 0;
        my $insb = 0;
        while($dna1 =~ m/(-+)/g) {
                $ins++;
                $insb+= length($&);
        }

        my $dels = 0;
        my $delb = 0;
	# dashes in dna2 indicate deletions in reference
	while ($dna2 =~ m/(-+)/g) {
                $dels++;
                $delb+= length($&);
	}

	my $mis  = 0;
        my $misb = 0;
        while($mat =~ m/\s+/g) {
                $mis++;
                $misb+= length($&);
        }

	$dna1 =~ s/-+//g;
        $dna2 =~ s/-+//g;

	$h->{ins}  = $ins;
	$h->{insb} = $insb;
	$h->{dels} = $dels;
	$h->{delb} = $delb;
	$h->{mis}  = $mis - $ins - $dels;
	$h->{misb} = $misb - $delb - $insb;
	$h->{cleandna1} = $dna1;
	$h->{cleandna2} = $dna2;

	#print "cleandna1 is $dna1\n";
	#print "cleandna2 is $dna2\n";


	return($h);
}
