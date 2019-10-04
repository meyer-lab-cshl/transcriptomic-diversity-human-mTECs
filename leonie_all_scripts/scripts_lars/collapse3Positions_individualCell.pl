#!/home/stroemic/software/miniconda2/envs/rnaseq/bin/perl

#purpose: define 3' sites OVER ONE SIGNLE CELL
#at the same time, produce an output that summarises counts & barcode counts per gene, then R doesn't need to do it later.

use strict;
use warnings;

my ($isoform_deffile, $out_file_isoforms, $out_file_genes) = @ARGV;

#part 1: read in the definitions of the isoforms. make a hash which for each gene stores the positions of called isoforms.
my ($line);
my %isoforms;
my @line;

open(ISOFORMS, "<$isoform_deffile");
while ($line = <ISOFORMS>) {
    @line = split(",", $line);
    $isoforms{$line[0]} = [] unless $isoforms{$line[0]};
    push(@{ $isoforms{$line[0]} }, $line[1] ); #store the isoform positions.
}

close(ISOFORMS);


#now copllapse to closest defined isoform. therefore, read in cell.
my ($add, $gene, $pos, $minpos, $dist, $mindist, $bccount, $readcount, $bc_antisense, $reads_antisense);
my $previous = "";
my @lines = ();


open(COLLAPSED, ">$out_file_isoforms");
open(GENES, ">$out_file_genes");
$line = <STDIN>;
print GENES "geneID,BarcodeCount,ReadCount,Barcode_antisense,Read_antisense\n";
print COLLAPSED $line;
while ($line = <STDIN>) {

    chomp $line;
    $add = [split(",",$line)]; #creates a pointer to the line.
    
    #even while reading in the file, change the position to the closest known isoform.
    #gene = $add->[0];
    #pos = $add->[1]; 
    
    #find the isoform with the closest distance to $gene/$pos
    $mindist = 1000; #changed
    $minpos = 0; #for each entry, first set the nearest position to 0
    foreach $pos (@{$isoforms{$add->[0]}}) {
        $dist = abs($add->[1] - $pos);
        if ($dist < $mindist) {
            $mindist = $dist;
            $minpos = $pos; #update the position if a nearby isoform has been found
        }
    }
    next unless $minpos; #skip the entry if no position has been found (if the value of $minpos is still 0)
    
    ###previous code; updated on 11.04.17 
    ##find the isoform with the closest distance to $gene/$pos
    #$mindist = 10000000;
    #foreach $pos (@{$isoforms{$add->[0]}}) {
    #    $dist = abs($add->[1] - $pos);
    #    if ($dist < $mindist) {
    #        $mindist = $dist;
    #        $minpos = $pos;
    #}
    #}

        $add->[1] = $minpos;
    
    if ($add->[0] ne $previous) { #if the gene name is not identical to the previous...
        collate(@lines); #collate all previous lines
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

collate(@lines); #collate all previous lines
print GENES "$previous,$bccount,$readcount,$bc_antisense,$reads_antisense\n" if $previous;
    
    
close(COLLAPSED);

#need to count # of barcodes at each position first!

sub collate {

    my @sorted;
    my ($current, $candidate, $merge, $dist, $k, $n, $p, $pval, $abstand);
    #STEP 3: Sort the 3' positions based on the number of barcodes found.
    @sorted = sort {$a->[2] <=> $b->[2]} @_;
    

    #STEP 5: Repeat the following three steps until the list of nonredundant sequences 
    #is empty
    #STEP 5A: Eliminate the bottom read of lowest abundance, S, from the list of nonredundant 
    #reads sorted by their frequencies in descending order.
    while ($current = shift @sorted) {

        $merge = 0; #assume the sequence will not be merged
        foreach $candidate (reverse @sorted) {
            $dist = abs($candidate->[1] - $current->[1]);

            if ($dist == 0) { 
                $merge = $candidate;
                last;
            }
        }
        
        #if a potential candidate has been found (same gene, similar 3' position)
        if ($merge) {
                #$merge->[1] = int(($merge->[1] * $merge->[2] + $current->[1] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
                #$merge->[5] = int(($merge->[5] * $merge->[2] + $current->[5] * $current->[2]) / ($merge->[2] + $current->[2]) + 0.5);
                $merge->[2] += $current->[2]; #add the barcode count to the parent
                $merge->[3] += $current->[3]; #add the read count to the parent
                
                } else {
           #call a 3' end if a sufficient number of barcodes has been identifed that support it
            print COLLAPSED join(",", @{$current}) . "\n"; #write the corrected line
        
  
            
        } 
        
    }

}


