#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);


my $usage = "converts a table of data to nexus format. States of 2 digits will be in {}\n\n./table2nexus.pl -in <FILE in table format with row names as languages> \n\n";

if($#ARGV < 0){ die $usage; }

my $input = "";
while( my $args = shift @ARGV){
    if($args =~ /^-in$/i){ $input = shift @ARGV; next; }
    die "Argument $args is invalid\n$usage";
}

my @species = ();
my @features = ();
my @v = ();
open(IN, $input) or die $!;

## read in the main data
my     $first = 1;
my $featcnt = 0;
my @data = ();
while( defined(my $ln = <IN>)){
    if($ln =~ /^\s*$/){ next; }
    chomp($ln);
    @v = split(/\s+/, $ln);
    push @species, shift @v;
    $featcnt = 0;
    my @strdata = ();

    foreach my $el (@v){
	$featcnt++;
	if($first){
	    push @features, $featcnt;
	}
	push @strdata, $el;
    }
    $first = 0;
    
    #print $datastr, "\n";
    ##push @data, $datastr;
    push @data, \@strdata;
}

for(my $i = 0; $i < @species; ++$i){
    print ">$species[$i]\n";
    my $curj = 0;
    
    for(my $j = 0; $j < @features; ++$j){
	print $data[$i][$j];
    }
    print "\n";
}
