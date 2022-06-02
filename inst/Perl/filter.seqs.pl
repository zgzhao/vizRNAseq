#!/usr/bin/env perl
# -*- coding:utf-8 -*-
# Description: 
# AUTHOR: ZG Zhao
# 2022-06-02 17:43:18
# Filter and prepare sequences from large file such as NCBI nt/nr for R

use utf8;
binmode STDIN, ':encoding(utf8)';
binmode STDOUT, ':encoding(utf8)';
use 5.010;
use File::Basename;

die "USAGE: $0 blast-fmt6-file fasta outfile\n" if $#ARGV<2;

## parse blastout: get unique matched names
open(iFH, "<$ARGV[0]") || die "Cannot open blast result.";
while(<iFH>){
  @para = /([^\s]+)/g;
  next if $#para < 2;
  $seq = $para[1];  ## matched seqname
  $blResult{$seq}++;
}
close(iFH);

@snames = keys %blResult;  ## array of matched names

## filter fasta file
open(iFH, "<$ARGV[1]");
open(oFH, ">$ARGV[2]");

$matched = 0;
while(<iFH>){
  next if(! $matched && /^[^>]/);
  if(/^>/) {
    $matched = 0;
    @name = /^>([^\s]+)/;
    foreach (@snames) {
      next if $matched;
      if($_ eq $name[0]) {
	$matched = 1;
      }
    }
  }
  ## print matched sequences only
  next if ! $matched;
  print oFH $_;
}
close(oFH);
close(iFH);

