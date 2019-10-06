#!~/anaconda3/envs/tss/bin/perl

# Author: Lars Velten, Hannah Meyer
# Purpose:
#   * define consensus sites within a sample
#   * produce an output that summarises counts & barcode counts per gene

###################
## subroutines ####
###################

sub usage {
    print "Unknown option: @_\n" if (@_);
    print "usage: Map sample transcription start sites (tss) to consensus tss";
    print "\t--insummary:\tSample-specific summary file of tss";
    print "\t--iniso:\tFile with consensus positions";
    print "\t--outiso:\tFile for sample-specific isoform output file";
    print "\t--outgene:\tFile for sample-specific gene output file";
    print "\t--h:\t";
    exit;
}


sub collate {
    my @sorted;
    my ($current, $candidate, $merge, $dist, $k, $n, $p, $pval, $abstand);
    # Sort positions based on the number of barcodes found.
    @sorted = sort {$a->[2] <=> $b->[2]} @_;

    # Repeat until list of nonredundant sequences is empty
    # Eliminate the bottom read of lowest abundance, S, from the list
    # of nonredundant reads sorted by their frequencies in descending order.
    while ($current = shift @sorted) {
        #assume the sequence will not be merged
        $merge = 0;
        foreach $candidate (reverse @sorted) {
            $dist = abs($candidate->[1] - $current->[1]);
            if ($dist == 0) {
                $merge = $candidate;
                last;
            }
        }
        # if a potential candidate has been found (same gene, similar position)
        if ($merge) {
            #add the barcode count to the parent
            $merge->[2] += $current->[2];
            #add the read count to the parent
            $merge->[3] += $current->[3];
        } else {
            # sufficient number of barcodes have been identifed that support it
            #write the corrected line
            print COLLAPSED join(",", @{$current}) . "\n";
        }
    }
}
##############################
## settings and variables ####
##############################

use strict;
use warnings;
use Getopt::Long;

my ($line, $add, $gene, $pos, $minpos, $dist, $mindist, $bccount, $readcount,
    $bc_antisense, $reads_antisense);
my %isoforms;
my @line;

# command line parameters
my ($iniso, $insummary, $outiso, $outgenes);

################
## analysis ####
################

if (@ARGV <  1 or ! GetOptions ("iniso=s" => \$iniso,
        "insummary=s" => \$insummary,
        "outiso=s" => \$outiso,
        "outgene=s" => \$outgene,
        "h=s" => \$h)) {
    usage()
}

GetOptions ("iniso=s" => \$iniso,
    "insummary=s" => \$insummary,
    "outiso=s" => \$outiso,
    "outgene=s" => \$outgene,
    "h=s" => \$h)

# command line parameters
print "Sample-specific summary file of tss: $insummary\n";
print "File with consensus positions: $iniso\n";
print "Path to file for sample-specific isoform output file: $outiso\n";
print "Path to file for sample-specific gene output file: $geneiso\n";

################
## analysis ####
################

# Read definitions of isoforms
# make a hash which for each gene stores the positions of called isoforms.
open(ISOFORMS, "<$iniso");
while ($line = <ISOFORMS>) {
    @line = split(",", $line);
    $isoforms{$line[0]} = [] unless $isoforms{$line[0]};
    # store the isoform positions.
    push(@{$isoforms{$line[0]}}, $line[1] );
}
close(ISOFORMS);


# now collapse to closest defined isoform. therefore, read in cell.
my $previous = "";
my @lines = ();

open(GENES, ">$outgenes");
print GENES "geneID,BarcodeCount,ReadCount,Barcode_antisense,Read_antisense\n";
open(COLLAPSED, ">$outiso");
print COLLAPSED "geneID,Position,BarcodeCount,ReadCount,PosFromAnno,Class\n";

open(SUMMARY, "<$insummary");
while ($line = <SUMMARY>) {
    chomp $line;
    # creates a pointer to the line.
    $add = [split(",",$line)];
    # find the isoform with the closest distance to $gene/$pos
    $mindist = 1000;
    # for each entry, first set the nearest position to 0
    $minpos = 0;
    foreach $pos (@{$isoforms{$add->[0]}}) {
        $dist = abs($add->[1] - $pos);
        if ($dist < $mindist) {
            $mindist = $dist;
            #update the position if a nearby isoform has been found
            $minpos = $pos;
        }
    }
    # skip the entry if no position has been found (if value $minpos still 0)
    next unless $minpos;
    $add->[1] = $minpos;
    #if the gene name is not identical to the previous...
    if ($add->[0] ne $previous) {
        #collate all previous lines
        collate(@lines);
        print GENES "$previous,$bccount,$readcount,$bc_antisense,$reads_antisense\n" if $previous;
        $bccount = 0;
        $readcount = 0;
        $bc_antisense = 0;
        $reads_antisense = 0;
        @lines = ();
        $previous = $add->[0];

    }
    $bccount += $add->[2] if ($add->[4] ne "antisense_gene") and ($add->[4] ne "intergenic") and ($add->[4] ne "antisense_upstream");
    $readcount += $add->[3] if ($add->[4] ne "antisense_gene") and ($add->[4] ne "intergenic") and ($add->[4] ne "antisense_upstream");
    $bc_antisense += $add->[2] if ($add->[4] eq "antisense_gene") or  ($add->[4] eq "antisense_upstream");
    $reads_antisense += $add->[3] if ($add->[4] eq "antisense_gene") or  ($add->[4] eq "antisense_upstream");
    push(@lines, $add);
}

# Collate all previous lines
collate(@lines);
print GENES "$previous,$bccount,$readcount,$bc_antisense,$reads_antisense\n" if $previous;
close(COLLAPSED);
