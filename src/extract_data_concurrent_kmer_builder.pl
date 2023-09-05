#!/bin/perl


# Extract system level

@all_logs = glob("./*.log");
%INPUT_SIZES;
%GENOME_NAME;
%CORE_COUNTS;
%ALGO_MODES;
%DATA;


foreach $f (@all_logs) {
  open(fptr, "<", $f) or die($!);
  $line = <fptr>;
  if (!($line =~/reading\s/)) {
    printf("...skipping file $f....\n");
    close(fptr);
    next;
  }

  print "processing $f\n";

  $gname; $algo_mode; $input_size; $core_count=1;
  $unique_kmers; $elapsed_time;
  if ($line =~ /([\w]+)\./) {
    $gname = $1;
  }

  while (<fptr>) {
    $line = $_;
    chomp($line);
    if ($line =~ /mode=([\w\-]+)/) {
      $algo  = $1;
      if ($algo eq "sort") {
        $algo = "Seq-Slide-Window-Sort";
      } elsif ($algo eq "hash_map") {
        $algo = "Seq-Slide-Window-Hash";
      } elsif ($algo eq "explicit") {
        $algo = "Concurrent-Slide-Window";
      } elsif ($algo eq "implicit-v0") {
        $algo = "Concurrent-Implicit-v0";
      } elsif ($algo eq "implicit-v1") {
        $algo = "Concurrent-Implicit-v1";
      } elsif ($algo eq "implicit") {
        $algo = "Concurrent-No-Slide-Window";
      }
    }

    if ($line =~ /tcount=([\d]+)$/) {
      $core_count = $1;
    }

    if ($line =~ /Reads\s*=\s*([\d]+)/) {
      $input_size = $1;
    }

    if ($line =~ /[Uu]nique [kK]mers=([\d]+)$/) {
      $unique_kmers = $1;
    }


    if ($line =~ /build elapsed time\]\]:\s+([\d\.]+)/) {
      $elapsed_time = $1;
    }

  }


  $key = $algo.":".$gname.":".$input_size.":".$core_count;
  $key_u = $key.":unique_kmers";
  $key_e = $key.":elapsed_time";


  $DATA{$key_u} = $unique_kmers;
  $DATA{$key_e} = $elapsed_time;

  $ALGO_MODES{$algo} = 1;
  $INPUT_SIZES{$input_size} = 1;
  $CORE_COUNTS{$core_count} = 1;
  $GENOME_NAME{$gname} = 1;

  close(fptr);
}

foreach $alg (@ALGO_MODES) {
  print "algo = $alg\n";
}

open(fptr, ">", "./full_system_results_07132023.csv") or die($!);

@system_cores = (1, 48);
@select_keys = ("unique_kmers", "elapsed_time");
#header #

print fptr "INPUT-SIZE, GENOME, ";

foreach $algo (reverse sort keys %ALGO_MODES) {
  foreach $select (@select_keys) {
    print fptr "$algo-$select,";
  }
}
print fptr "\n";


foreach $input_size (sort {$a <=> $b } (keys %INPUT_SIZES)) {
  foreach $gname (keys %GENOME_NAME) {
    $row_start = 0;
    foreach $core_count (@system_cores) {
      foreach $algo (reverse sort keys %ALGO_MODES) {
        foreach $select (@select_keys) {
          $key = $algo.":".$gname.":".$input_size.":".$core_count;
          $key_s = $key.":".$select;
          if (exists $DATA{$key_s}) {
            if ($row_start == 0) {
              printf fptr "$input_size, $gname, ";
              $row_start = 1;
            }
            $v = $DATA{$key_s};
            printf fptr "$v, ";
          }
        }
      }
    }
    if ($row_start == 1) { print fptr "\n"; }
  }
}
close(fptr);


open(fptr, ">", "./thread_scaling_concurrent_algos_results_07132023.csv") 
  or die("$!");

@study_modes = ("Concurrent-Implicit-v0", "Concurrent-Implicit-v1");

$max_size =0;
foreach $key (keys %INPUT_SIZES) {
  if ($key > $max_size) {
    $max_size = $key;
  }
}
printf fptr "CORE-COUNT, INPUT-SIZE";
foreach $algo (@study_modes) {
  printf fptr ",$algo-elapsed_time(sec)";
}
printf fptr "\n";

$select_gname = "chimp";
foreach $core_count (sort {$a <=> $b} (keys %CORE_COUNTS)) {
  $row_start = 0;
  foreach $algo (@study_modes) {
    $key = $algo.":".$select_gname.":".$max_size.":".$core_count;
    $key = $key.":elapsed_time";
    if (exists $DATA{$key}) {
      if ($row_start == 0) {
        printf fptr "$core_count, $max_size, ";
        $row_start = 1;
      }
      $v = $DATA{$key};
      printf fptr "$v,";
    }
  }
  printf fptr "\n";
}
close(fptr);




