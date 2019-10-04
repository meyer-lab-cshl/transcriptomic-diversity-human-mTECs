#!/home/stroemic/software/miniconda2/envs/rnaseq/bin/perl

#purpose: define 3' sites OVER ALL CELLS.
#empirically determined variation of 3' ends (R-script): 11 -collapse all ends that are in a 11bp window surrounding the max. peak

#the window is not strictly 11 bp. Because it will merge read 1 to read 2 that is within an 11bp wndow. Then it will merge read 2 to sth that is in a 11bp window. etc. 
#Starting with the read with the smallest number of barcodes. It will always merge to the read with the largest number of barcodes within that window.

use strict;
use warnings;

#but first, read in all count tables
print "now in perl script\n";

my ($count_dir, $outfile_collapsed, $outfile_summarised, $max_distance, $min_no_bc) = @ARGV;

print "$count_dir\n";
print "$outfile_collapsed\n";
print "$outfile_summarised\n";
print "$max_distance\n";
print "$min_no_bc\n";

my @files = glob($count_dir . "*.summary.counts");

glob($count_dir . "*.summary.counts");
#print join(",", @files);

my @lines;
my ($line, $file, $add, $cell, $header);

foreach $file (@files) {
    print "Now working on $file \n";
    open(READER, "<$file"); #for debugging: read in only ENSMUSG00000024067
    #open(READER, "grep ENSMUSG00000024067 $file|");
    die "could not identify cell name" unless ($file =~ /.+\/(.+)\.summary\.counts/);
    $cell = $1; 
    $header = <READER>;
    while ($line = <READER>) {
        chomp $line;
        
        $add = [split(",",$line) ]; #creates a pointer to the line.
        push(@lines, $add);
    }
    close(READER)
}

#all count tables are now stored in @lines

@lines = sort {$a->[0] cmp $b->[0] } @lines;
my @current_gene = ();


open(my $collapseFH, ">$outfile_collapsed");
open(my $summariseFH, ">$outfile_summarised");
my $previous = "";
my ($current, $summarised, $final);
while ($current = shift @lines) {
    unless ($current->[0] eq $previous) {
        #before calling collate, I need to count the barcodes per position! without merging nearby positions!
        #collate(@current_gene, 1,0);
        #print $current->[0] . "\n";
        if ($#current_gene > 50000) {warn("$previous has too many entries, omitting\n")} else { 
        $summarised = collate(\@current_gene, 0,$summariseFH); #sumarise all previous lines - just over all cells count barcodes 
        $final = collate($summarised, $max_distance,$collapseFH);
        }
        #the result is now in %result. print it.
        @current_gene = ();
        $previous = $current->[0];
    }
    push(@current_gene, $current);
}

$summarised = collate(\@current_gene, 0,$summariseFH); #sumarise all previous lines - just over all cells count barcodes 
$final = collate($summarised, $max_distance,$collapseFH);

close($collapseFH);
close($summariseFH);


#need to count # of barcodes at each position first!
#for debugging, report all combinations of cells / pos for gene ENSMUSG00000024816

sub collate {
    my $md = $_[1]; #allows passing an option: what is the max distance to merge? run with 0 
    my $write_out = $_[2]; #2nd option: should output be written?
    my @sorted;
    my @out; #return value: the collapsed array
    my ($current, $candidate, $merge, $dist, $k, $n, $p, $pval, $abstand, $rem);
    #STEP 3: Sort the 3' positions based on the number of barcodes found.
    @sorted = sort {$a->[2] <=> $b->[2]} @{$_[0]};
    #foreach $candidate (reverse @sorted) {
    #print $candidate->[1] . "\t" . $candidate->[2] . "\n"; #debugging only.
    #}
    #print $current->[1] . "\t" . $current->[2]; #debugging only.

    #STEP 5: Repeat the following three steps until the list of nonredundant sequences 
    #is empty
    #STEP 5A: Eliminate the bottom read of lowest abundance, S, from the list of nonredundant 
    #reads sorted by their frequencies in descending order.
    while ($current = shift @sorted) {
    #print $current->[1] . "\t" . $current->[2]; #debugging only.
        $merge = 0; #assume the sequence will not be merged
        foreach $candidate (reverse @sorted) {
            
            $dist = abs($candidate->[1] - $current->[1]);

            if ($dist <= $md) { 
                $merge = $candidate;
                last;
            }
        }
        if ($current->[0] eq "ENSMUSG00000024816") {
            $rem = $merge ? $merge->[1] : 0;
            print "$current->[1], $rem \n";
    }
        #if a potential candidate has been found (same gene, similar 3' position)
        if ($merge) {
                $merge->[1] = int(($merge->[1] * $merge->[2] + $current->[1] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
                $merge->[5] = int(($merge->[5] * $merge->[2] + $current->[5] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
                $merge->[2] += $current->[2]; #add the barcode count to the parent
                $merge->[3] += $current->[3]; #add the read count to the parent
                
                } elsif ($current->[2] >= $min_no_bc) {
    #call a 3' end if a sufficient number of barcodes has been identifed that support it
            print $write_out join(",", @{$current}) . "\n" if $write_out; #write the corrected lie
            push(@out, $current); #store the final line
            
        } 
        
    }

 return \@out;
}

