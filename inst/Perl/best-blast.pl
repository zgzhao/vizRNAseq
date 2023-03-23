#!/usr/bin/env perl
# -*- coding:utf-8 -*-
# Description: get best blast accessions
# AUTHOR: ZG Zhao
# 2022-07-06 18:32:31

use utf8;
use strict;
binmode STDIN, ':encoding(utf8)';
binmode STDOUT, ':encoding(utf8)';
use v5.10.1;
use List::Util qw(min max);

sub uniq {
  my %temp;
  return grep { !$temp{$_}++ } @_;
}

sub nsort {
  sort {$a <=> $b} @_;
}

sub nrange {
  my @ans;
  while(@_) {
    my $ss = shift @_;
    my $ee = shift @_;
    $ss = $ss/100;
    $ee = $ee/100;
    @ans = uniq(@ans, ($ss..$ee));
  }
  my $len = @ans;
  return $len;
}
## qid, sid, p_indentical, n_align, n_mismatch, gap_open
## qstart, qend, sstart, send, eval, score

open iFH, "<$ARGV[0]";
my %ranges;
my %scores;
while(<iFH>) {
  chomp;
  my @info = /([^\s]+)/g;
  my $len = @info;
  next if $len != 12;
  my $xid  = $info[0]. " " .$info[1];
  my $rnx = min($info[6], $info[7])." ".max($info[6], $info[7]);
  if(exists $ranges{$xid}) {
    @ranges{$xid} .= " $rnx";
    $scores{$xid} += $info[11];
  } else {
    @ranges{$xid} = $rnx;
    $scores{$xid} = $info[11];
  }
}
close iFH;

my %nmatch;
foreach my $xid (keys %ranges) {
  $nmatch{$xid} = nrange(split(" ", %ranges{$xid}));
}
%ranges = ();

my %bestid;
my %bestmm;
my %bestss;
foreach my $xid (keys %nmatch) {
  my ($qid, $sid) = split(" ", $xid);
  my $mm = $nmatch{$xid};
  my $bb = $scores{$xid};
  next if($mm < $bestmm{$qid});
  if($mm > $bestmm{$qid} || $bb > $bestss{$qid}) {
    $bestid{$qid} = $sid;
    $bestmm{$qid} = $mm;
    $bestss{$qid} = $bb;
  }
}

open oFH, ">$ARGV[1]";
foreach my $qid (keys %bestid) {
  print oFH "$qid\t$bestid{$qid}\t$bestmm{$qid}\t$bestss{$qid}\n"
}
close oFH;

