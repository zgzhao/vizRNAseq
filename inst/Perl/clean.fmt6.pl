#!/usr/bin/env perl
# -*- coding:utf-8 -*-
# Description: 
# AUTHOR: ZG Zhao
# 2022-06-02 17:43:18
# Filter and prepare sequences from large file such as NCBI nt/nr for R

## TODO: to be implemented.

use utf8;
binmode STDIN, ':encoding(utf8)';
binmode STDOUT, ':encoding(utf8)';
use 5.010;
use File::Basename;

die "USAGE: $0 blast-fmt6-file outfile\n" if $#ARGV<1;

open(iFH, "<$ARGV[0]") || die "Cannot open blast result.";
open(oFH, ">$ARGV[1]");

$lid = "";
while(<iFH>){
  @minfo = /([^\s]+)/g;
  next if $#minfo < 10;
  $qid   = $minfo[0];

  if ($qid ne $lid) {
    %sinfo = {};
  }
  ## 按目标序列整理匹配区域
  $sid   = $minfo[1];
  $sndx1 = $minfo[8];
  $sndx2 = $minfo[9];

  ##
  next if ! $matched;
  ## 输出最佳匹配结果
  print oFH $_;
}

close(iFH);
close(oFH);

