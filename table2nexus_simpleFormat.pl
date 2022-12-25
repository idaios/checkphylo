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

print "#NEXUS\n\n";

print "BEGIN DATA;\n";
print "\tDIMENSIONS NTAX=", scalar @species, " NCHAR=", scalar @features, ";\n";
print "\tFORMAT DATATYPE=Standard SYMBOLS= \"01\" MISSING=? GAP= -;\n";
print "\tCHARSTATELABELS\n";

# print "\tTAXLABELS\n";

# foreach my $s (@species){
#     print "\t\t$s\n";
# }
# print "\t\t;\n";
# print "END;\n\n\n";

# print "BEGIN CHARACTERS;\n";
# print "[Section for the data follows]\n\n\n";

# print "\tTITLE Appendix_1;\n";
# print "\tLink TAXA = combined;\n";
# print "\tDIMENSIONS NCHAR=",$featcnt, ";\n";
# print "\tFORMAT DATATYPE=Standard SYMBOLS= \"01\" MISSING=? GAP= -;\n";
# print "\tCHARSTATELABELS\n";

# my @uniqueFeatures = ();

# my @newmapfeatures = ();
# my @classes = ();
# my $divClass = 5;
# my $maxClass = 9;
# push @classes, "absent";

# for( my $i=0; $i < $maxClass-1; ++$i){
#     push @classes, "bet_".$i*$divClass."-".($i+1)*$divClass;
# }
# push @classes, "greaterThan_".($divClass)*($maxClass-1);

# for (my $i = 0; $i < @features; ++$i){
#     my @tmpv = uniq @{$featurestable[$i]};
#     my @v = ();
#     for(my $j = 0; $j < @{$featurestable[$i]}; $j++){
# 	my $cl = int($featurestable[$i]/5);
#     }
# }



for( my $i = 0; $i < $featcnt; ++$i){
    print "\t\t",$i+1," ";
    print "'",$features[$i],"'"; ##,","; ##"' / ";
    ##my @tmpv = uniq @{$featurestable[$i]};
    ##print join('. ',  uniq @{$featurestable[$i]});
##    print join(" ", @classes);
    # my @tmpv = uniq @{$featurestable[$i]};
    # my %tmphash = ();
    # for(my $j = 0; $j < @tmpv; ++$j){
    # 	$tmphash{$tmpv[$j]} = $j;
    # }
    # push @uniqueFeatures, \%tmphash;
    print ",\n";
}
print "\t\t;\n";

print "MATRIX\n";
for(my $i = 0; $i < @species; ++$i){
    print $species[$i], "\t";
    my $curj = 0;
    
    for(my $j = 0; $j < @features; ++$j){
	print $data[$i][$j];
    }
	# }else{
    # 	for(my $j = 0; $j < @stabfeatures; ++$j){
    # 	    for(my $jj = 0; $jj < @{$stabilityFeatures[$j]}; ++$jj){
    # 		my $index = $stabilityFeatures[$j][$jj];
    # 		if($jj == 0){print "0";}
    # 		print $data[$i][$index];
    # 	    }
    # 	}
    # }
    print "\n";
}

print ";\nEND;\n";

# print "begin assumptions;\n";
# if($stabilityClasses == 0){
#     for(my $i = 0; $i < @ordlan; ++$i){
# 	print "\tcharset ", $ordlan[$i], " = ", join(" ", @{$langFeatures[$i]}), ";\n";
#     }
# }else{
#     my $curj = 1;
#     for(my $i = 0; $i < @stabfeatures; ++$i){
# 	print "\tcharset stab".$i." = ".$curj++;
# 	for(my $j = 0; $j < @{$stabilityFeatures[$i]}; ++$j){
# 	    print " ",$curj++;
# 	}
# 	print ";\n";
#     }
# }
print "\n";
print "end;\n";
    
    
    
    
