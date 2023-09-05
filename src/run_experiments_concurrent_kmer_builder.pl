#!/bin/perl

@files = glob("./*.fasta");
@core_count = (2, 4, 8, 10, 16, 18, 24, 32, 48);
foreach $f (@files) {

  $file_prefix;
  if ($f =~/.?\/?(.*?)\.fasta/) {
    $file_prefix = $1;
  } else { next; }

  print "$file_prefix\n";

  $run = "/usr/bin/time -o default.$file_prefix.time ./default-kmer-table-builder $f >& default.$file_prefix.log";
  printf("running: $run\n");
  system("$run");

  $run = "/usr/bin/time -o default.sort.$file_prefix.time ./default-kmer-table-builder $f --use_sort >& default.sort.$file_prefix.log";
  printf("running: $run\n");
  system("$run");

  foreach $core (@core_count) {
    $run1 = "/usr/bin/time -o concurrent-e.$file_prefix.$core.time ./concurrent-kmer-table-builder $f $core --use_explicit >& concurrent-e.$file_prefix.$core.log";

    $run2 = "/usr/bin/time -o concurrent-i.$file_prefix.$core.time ./concurrent-kmer-table-builder $f $core  >& concurrent-i.$file_prefix.$core.log";

    printf("$run1\n");
    printf("running....");
    system("$run1");

    printf("$run2\n");

    printf("running....");
    system("$run2");
 }

}

