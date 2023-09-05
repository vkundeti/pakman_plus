#!/bin/perl

@data_files = glob("./*.out");
$NEW_ALGO_NAME = "concurrent-kmer.v0";

#key: readcount-nodecount-algo-total-elapsed , "kmer-count", "start-modes", "end-mnodes", "uniq-kmers"


%DATA_HASH;
%READ_COUNTS;
%NODE_COUNTS;
%ALGOS;

foreach $f (@data_files) {
  open(FPTR, "<", $f) or die("$!");

  $mode = "out-of-box";

  if ($f =~ /concurrent-kmer\.v0/)  {
    $mode = "new-concurrent-algo-v0";
  } elsif ($f =~ /concurrent-kmer\.v1/) {
    $mode = "new-concurrent-algo-v1";
  }


# read-count
# node-count
# total-elapsed
# kmer-count
# start-mnodes
# end-mnodes
# uniq-kmers 

  $read_count=""; $node_count=""; 
  $total_elapsed=""; $kmer_count_time="";
  $start_mnodes=""; $end_mnodes="";
  $uniq_kmers="";

  while (<FPTR>) {
    $line = $_;

    if ($line =~ /Total number of reads:\s+([\d]+)/) {
      $read_count = $1;
    }

    if ($line =~ /,\s+Number of Processes:\s+([\d]+)/) {
      $node_count = $1;
    }

    if ($line =~ 
          /Average time for performing k-mer counting[^:]+:\s+([\d\.]+)/) {
      $kmer_count_time = $1;
    }

    if ($line =~ /Total distinct k-mer entries[^:]+:\s+([\d]+)/) {
      $uniq_kmers = $1;
    }

    if ($line =~ /Itr: 1, Total number of macro_nodes across[^:]+: ([\d]+)/) {
      $start_mnodes = $1;
    }

    if ($line =~ /Itr: [\d]+, Total number of macro_nodes across[^:]+: ([\d]+)/) {
      $end_mnodes = $1;
    }

    if ($line =~ /elapsed-time-no-io\]:\s+([\d\.]+)/) {
      $total_elapsed = $1;
    }
  }

  $key = $mode.":".$read_count.":".$node_count;
#  print "$mode, $f, $read_count, $node_count, $kmer_count_time, ";
#  print "$start_mnodes, $end_mnodes, $total_elapsed\n";

  $ALGOS{$mode} = 1;
  $READ_COUNTS{$read_count} = 1;
  $NODE_COUNTS{$node_count} = 1;
 
  if (not ($kmer_count_time eq "")) {
    $DATA_HASH{$key.":kmer_count_time"} = $kmer_count_time;
  }

  if (not ($start_mnodes eq "" ) ) {
    $DATA_HASH{$key.":start_mnodes"} = $start_mnodes;
  }

  if (not ($end_mnodes eq "" )) {
    $DATA_HASH{$key.":end_mnodes"} = $end_mnodes;
  }

  if (not ($total_elapsed eq "")) {
    $DATA_HASH{$key.":total_elapsed"} = $total_elapsed;
  }

}

#@modes_of_interest = ("out-of-box", "new-concurrent-algo-v0", "new-concurrent-algo-v1");

@modes_of_interest = ("out-of-box",  "new-concurrent-algo-v1");
@select_list = 
  ("start_mnodes", "end_mnodes", "kmer_count_time", "total_elapsed");


##########################header###############################
{
{
    
    print "read_count, node_count,";
    foreach $mode (@modes_of_interest) {
      foreach $select (@select_list) {
        print "$mode-$select,";
      }
    }
    print "\n";
} # foreach node #
} # foreach read #
###############################################################

##########################data#################################
foreach $read_count (sort {$a <=> $b } (keys %READ_COUNTS)) {
  foreach $node_count (sort {$a <=> $b} (keys %NODE_COUNTS)) {
    print "$read_count, $node_count,";

    foreach $mode (@modes_of_interest) {
      $key = $mode.":".$read_count.":".$node_count;
      foreach $select (@select_list) {
        $skey = $key.":".$select;
        $value = "";
        if (exists $DATA_HASH{$skey}) {
          $value = $DATA_HASH{$skey};
        }
        print "$value, ";
      }
    }
    print "\n";
  } # foreach node #
} # foreach read #
