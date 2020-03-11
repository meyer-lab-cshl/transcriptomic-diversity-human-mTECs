#!~/anaconda3/envs/tss/bin/perl

# Author: Lars Velten, Hannah Meyer
# Purpose: define consensus sites over all samples.
# collapse all counts in a $maxdist bp window surrounding the max. peak

# the window is not strictly $maxdist bp:
# -> it will merge read 1 to read 2 that is within $maxdist window
# -> it will merge read 2 to sth that is in $maxdist window. etc.
# Starting with the read with the smallest number of barcodes:
# -> it will always merge to the read with the largest number of barcodes within
# that window.

###############
# subroutines #
###############

sub usage {
    print "Unknown option: @_\n" if (@_);
    print "usage: Collapse tss sites across experiment into consensus tss";
    print "\t--dir:\tDirectory(ies) with count files ending in --suffix; if multiple directories specified, separate by commas, no space
           --collapsed:\t
           --summarised:\t
           --maxdist:\tmaximum peak distance
           --minbc:\tminimum number of counts per position
           --suffix:\tsuffix of count files to be read";
    exit;
}

sub collate {
    # allows passing an option: what is the max distance to merge? run with 0
    my $md = $_[1];
    # 2nd option: should output be written?
    my $write_out = $_[2];
    my $bcfilter = $_[3];
    my @sorted;
    # return value: the collapsed array
    my @out;
    my ($current, $candidate, $merge, $dist, $k, $n, $p, $pval, $abstand, $rem);
    # Sort the positions based on the number of barcodes found.
    @sorted = sort {$a->[2] <=> $b->[2]} @{$_[0]};
    # Repeat until list of nonredundant sequences is empty
    # Eliminate the bottom read of lowest abundance, S, from the list
    # of nonredundant reads sorted by their frequencies in descending order.
    while ($current = shift @sorted) {
        # assume the sequence will not be merged
        $merge = 0;
        foreach $candidate (reverse @sorted) {
            $dist = abs($candidate->[1] - $current->[1]);
            if ($dist <= $md) {
                $merge = $candidate;
                last;
            }
        }
        #if a potential candidate has been found (same gene, similar position)
        if ($merge) {
            $merge->[1] = int(($merge->[1] * $merge->[2] + $current->[1] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
            $merge->[5] = int(($merge->[5] * $merge->[2] + $current->[5] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
            # add the barcode count to the parent
            $merge->[2] += $current->[2];
            # add the read count to the parent
            $merge->[3] += $current->[3];
        } elsif ($current->[2] >= $bcfilter) {
            # sufficient number of barcodes have been identifed that support it
            # write the corrected line
            print $write_out join(",", @{$current}) . "\n" if $write_out;
            push(@out, $current); #store the final line
        }
    }
 return \@out;
}

##############################
## settings and variables ####
##############################

use strict;
use warnings;
use Getopt::Long;

my @lines;
my ($line, $file, $add, $cell, $header);

# command line arguments
my ($dir, $collapsed, $summarised);

# command line arguments defaults
my $suffix = "summary.counts";
my $maxdist = 12;
my $minbc = 10;

################
## analysis ####
################

if (@ARGV <  1 or ! GetOptions ("dir=s" => \$dir,
        "collapsed=s" => \$collapsed,
        "summarised=s" => \$summarised,
        "suffix=s" => \$suffix,
        "maxdist=i" => \$maxdist,
        "minbc=i" => \$minbc)) {
    usage();
}

GetOptions (
    "dir=s" => \$dir,
    "collapsed=s" => \$collapsed,
    "summarised=s" => \$summarised,
    "suffix=s" => \$suffix,
    "maxdist=i" => \$maxdist,
    "minbc=i" => \$minbc);

# command line parameters
print "Directory with count files: $dir\n";
print "$collapsed\n";
print "$summarised\n";
print "Maximum peak distance: $maxdist\n";
print "Minimum counts per position: $minbc\n";
print "Suffix of counts files: $suffix\n";

# collect all count files from comma-separated directory list
my @dir = split(",", $dir);
my @files = ();
foreach my $d (@dir) {
    push(@files, glob($d . "/" . "*". $suffix));
}

foreach $file (@files) {
    print "Now working on $file \n";
    open(READER, "<$file");
    die "could not identify cell name" unless ($file =~ /.+\/(.+)$suffix/);
    $cell = $1;
    $header = <READER>;
    while ($line = <READER>) {
        chomp $line;
        # creates a pointer to the line.
        $add = [split(",",$line)];
        push(@lines, $add);
    }
    close(READER)
}

# all count tables are now stored in @lines
@lines = sort {$a->[0] cmp $b->[0] } @lines;
my @current_gene = ();

open(my $collapseFH, ">$collapsed");
open(my $summariseFH, ">$summarised");
my $previous = "";
my ($current, $summary, $final);
while ($current = shift @lines) {
    unless ($current->[0] eq $previous) {
        # before calling collate, count barcodes per position without merging
        # nearby positions!
        if ($#current_gene > 50000) {
            warn("$previous has too many entries, omitting\n")
        } else {
            # sumarise all previous lines - just over all cells count barcodes
            $summary = collate(\@current_gene, 0, $summariseFH, $minbc);
            $final = collate($summary, $maxdist, $collapseFH, $minbc);
        }
        # the result is now in %result. print it.
        @current_gene = ();
        $previous = $current->[0];
    }
    push(@current_gene, $current);
}

# sumarise all previous lines - just over all cells count barcodes
$summary = collate(\@current_gene, 0, $summariseFH, $minbc);
$final = collate($summary, $maxdist, $collapseFH, $minbc);

close($collapseFH);
close($summariseFH);
