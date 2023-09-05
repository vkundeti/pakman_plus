#!/bin/perl

@input_files = glob("./*.fasta");

foreach $f (@input_files) {
  print "$f\n";
}

$min_nodes=4;
$max_nodes=224;
$MAX_MODES = 3;


$run_id = 1;
foreach $f (@input_files) {
  $fprefix;
  if ($f =~ /\.\/(.*?)\.fasta/) {
    $fprefix = $1;
  } else {
    print "skipping $f\n";
    next;
  }
  $nodes = $min_nodes;
  while ($nodes <= $max_nodes) {
    for ($i=0; $i<$MAX_MODES; $i++) {
      $mode = "default";
      $mode_args = "";
      if ($i == 1) {
        $mode = "concurrent-kmer.v0";
        $mode_args = " -i -p 0 ";
      } elsif ($i == 2) {
        $mode = "concurrent-kmer.v1";
        $mode_args = " -i -p 1 ";
      }

      $cmd = "srun -p regular --nodes=$nodes --ntasks-per-node=1 --exclusive --output=pakman.$fprefix.$mode.$nodes-%J.out --error=pakman.$fprefix.$mode.$nodes-%J.err ./pakman-mpi -f $f -r 100 -c 100 -b 1000000000 $mode_args >& run.$run_id.log &";

      system($cmd);
      print "\n $cmd\n";
      $run_id = $run_id+1;
    }

    if ($nodes == $max_nodes) { last; }

    $nodes = $nodes*2;
    if ($nodes > $max_nodes) {
      $nodes = $max_nodes;
    }
  }
}
