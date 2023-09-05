#!/bin/perl

@log_files = glob("*.1Billion*default*.out");

$kmer_count_avg_percent=0;
$build_graph_avg_percent=0;
$itr_compact_avg_percent=0;
$gather_walk_avg_percent=0;
%DATA;
$sample_count = 0;
$read_count = 0;
$kmers = 0;
foreach $f (@log_files) {
  open(FPTR, "<", $f) or die($!);
  $node_count = 0;
  $kmers = 0;
  $kmer_count_time=0;
  $build_graph = 0;
  $itr_compact = 0;
  $itr_count=0;
  $total = 0;
  while (<FPTR>) {
    $line = $_;
    chomp($line);

    if ($line =~ /Total number of reads: ([\d]+)/) {
      $read_count = $1;
    }
    
    if ($line =~ /Average time for performing k-mer counting across all procs \(secs\): ([\d\.]+)/) {
      $kmer_count_time = $1;
    }
    if ($line =~ /Number of Processes: ([\d]+)/){
      $node_count = $1;
    }
    if ($line =~ /Total distinct k-mer entries across all proc\'s: ([\d]+)/) {
      $kmers = $1;
    }
    if ($line =~ /Average time for MN node construction across all procs \(secs\): ([\d\.]+)/) {
      $build_graph = $1;
    }

    if ($line =~ /Itr: ([\d]+), Total time for this iteration: ([\d\.]+)/ ) {
      $itr_compact += $2;
      $itr_count = $1; 
    }

    if ($line =~ /elapsed-time-no-io\]: ([\d\.]+)/) {
      $total_time = $1;
    }
  }
  close(FPTR);

  $all_except = $build_graph + $kmer_count_time + $itr_compact;
  $gather_walk = ($total_time - $all_except);

  $kmer_count_percent = $kmer_count_time/$total_time;
  $kmer_count_avg_percent += $kmer_count_percent;

  $build_graph_percent = $build_graph/$total_time;
  $build_graph_avg_percent += $build_graph_percent;

  $itr_compact_percent = $itr_compact/$total_time;
  $itr_compact_avg_percent += $itr_compact_percent;

  $gather_walk_percent = $gather_walk/$total_time;
  $gather_walk_avg_percent += $gather_walk_percent; 

  $build_graph_print = sprintf("%.3f", $build_graph);
  $itr_compact_print = sprintf("%.3f", $itr_compact);
  $kmer_count_print = sprintf("%.3f", $kmer_count_time);
  $gather_walk_print = sprintf("%.3f", $gather_walk);

  $build_graph_percent_print = sprintf("%.3f", ($build_graph/$total_time)*100);
  $itr_compact_percent_print = sprintf("%.3f", ($itr_compact/$total_time)*100);
  $kmer_count_percent_print = sprintf("%.3f", ($kmer_count_time/$total_time)*100);
  $gather_walk_percent_print = sprintf("%.3f", ($gather_walk/$total_time)*100);
  ++$sample_count;


# node, k-mers, $itr_count, kmer_count_time, $build_graph, $itr_compact, $gather_walk, $total, $itr_count #
  $DATA{$node_count} = " $itr_count & $kmer_count_print & $kmer_count_percent_print\\% & $build_graph_print & $build_graph_percent_print\\%  & $itr_compact_print & $itr_compact_percent_print\\% & $gather_walk_print & $gather_walk_percent_print\\%  & $total_time ";


}


# print table
print "\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|}\n";
print "\\hline\n";
print "\\multicolumn{11}{|c|}{Read Count=$read_count  Distinct k-mers=$kmers}\\\\\n";
print "\\hline\n";
print " \\multicolumn{2}{|c|}{} & \\multicolumn{2}{|c|}{Gen.k-mers}& \\multicolumn{2}{|c|}{GraphBuild} & \\multicolumn{2}{|c|}{Itr.Compact} & \\multicolumn{2}{|c|}{Gather-Walk} & End-to-End\\\\\n";
print "\\hline\n";
print " Nodes & Itr. Count & Time(sec) & Contrib. & Time(sec)& Contrib. & Time(sec) & Contrib. & Time(sec) & Contrib & Total-Time(sec) \\\\\n";
foreach $k ( sort { $a <=> $b } keys %DATA) {
  $v = $DATA{$k};
  print "\\hline\n";
  print "$k & $v \\\\ \n";
  
}

$kmer_count_avg_percent = $kmer_count_avg_percent/$sample_count;
$kmer_count_avg_print = sprintf("%0.3f", $kmer_count_avg_percent*100);

$build_graph_avg_percent = $build_graph_avg_percent/$sample_count;
$build_graph_avg_print = sprintf("%0.3f" , $build_graph_avg_percent*100);

$itr_compact_avg_percent = $itr_compact_avg_percent/$sample_count;
$itr_compact_avg_print = sprintf("%0.3f", $itr_compact_avg_percent*100);

$gather_walk_avg_percent = $gather_walk_avg_percent/$sample_count;
$gather_walk_avg_print = sprintf("%0.3f", $gather_walk_avg_percent*100);

print "\\hline\n";
print "\\multicolumn{11}{|c|}{Normalized Avg. Contribution of Component to Runtime}\\\\\\hline \n";
print " \\multicolumn{2}{|c|}{} & \\multicolumn{2}{|c|}{Gen.k-mers}& \\multicolumn{2}{|c|}{GraphBuild} & \\multicolumn{2}{|c|}{Itr.Compact} & \\multicolumn{2}{|c|}{Gather-Walk} & End-to-End\\\\\n";
print "\\hline\n";
print "\\multicolumn{2}{|c|}{} & \\multicolumn{2}{|c|} {$kmer_count_avg_print\\%} & \\multicolumn{2}{|c|}{$build_graph_avg_print\\%} & \\multicolumn{2}{|c|}{$itr_compact_avg_print\\%} & \\multicolumn{2}{|c|}{$gather_walk_avg_print\\%} & 100\\% \\\\\n";
print "\\hline\n";
print "\\end{tabular}\n";



