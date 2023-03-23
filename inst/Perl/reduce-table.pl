#!/usr/bin/env perl
# -*- coding:utf-8 -*-
# Description: BLAST2GO step2 - filter GeneID
# AUTHOR: ZG Zhao
# 2022-07-03 18:20:20

use utf8;
binmode STDIN, ':encoding(utf8)';
binmode STDOUT, ':encoding(utf8)';
use 5.010;
use File::Basename;

die "USAGE: $0 <file-ids> <file-maps> outfile\n" if $#ARGV<2;
die "$ARGV[0] not exists!" if(! -e $ARGV[0]);
die "$ARGV[1] not exists!" if(! -e $ARGV[1]);

open(iFH, "<$ARGV[0]");
while(<iFH>) {
  chomp;
  s /\.[0-9]+$//;
  $xids{$_} = 1;
}
close(iFH);

open(iFH, "<$ARGV[1]");
open(oFH, ">$ARGV[2]");

## columns: nr_protein_id, ncbi_gene_id
while(<iFH>){
  @minfo = /([^\s]+)/g;
  $acc = $minfo[0];
  $acc =~ s /\.[0-9]+$//;
  next if ! $xids{$acc};
  print oFH "$acc\t$minfo[1]\n";
}

close(iFH);
close(oFH);
