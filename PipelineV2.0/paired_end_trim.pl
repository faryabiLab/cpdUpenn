#!/usr/bin/perl -w

#
# a simple program for post-processing paired-end reads after btrim trimming
#
# Yong Kong
# Yale University 
#

use strict;
 
my $f1 = shift;  # summary file 1
my $f2 = shift;  # summary file 2
my $s1 = shift;  # trim output file 1
my $s2 = shift;  # trim output file 2

my $check = shift || 0;

die "usage: $0 <summary file 1> <summary file 2> <trim output file 1> <trim output file 2>" unless ($f1 && $f2 && $s1 && $s2);

my ($fh1, $fh2);
my ($sh1, $sh2);
my ($oh1, $oh2);
my ($ph1, $ph2);

open($fh1, "<$f1") || die "cannot open $f1";
open($fh2, "<$f2") || die "cannot open $f2";

open($sh1, "<$s1") || die "cannot open $s1";
open($sh2, "<$s2") || die "cannot open $s2";

my $o1 = "${s1}.pe";  # reads in both ends passed
my $o2 = "${s2}.pe";

my $p1 = "${s1}.se";  # read in trim file1 passed, but the pair in file2 failed
my $p2 = "${s2}.se";  # read in trim file2 passed, but the pair in file1 failed

open($oh1, ">$o1") || die "cannot open $o1";
open($oh2, ">$o2") || die "cannot open $o2";

open($ph1, ">$p1") || die "cannot open $p1";
open($ph2, ">$p2") || die "cannot open $p2";

while (my $l1 = <$fh1>, my $l2 = <$fh2>) {
    chomp($l1);
    my @t1 = split(/\t/, $l1);
    my $n1 = $t1[0];
    my $s1 = $t1[1];

    chomp($l2);
    my @t2 = split(/\t/, $l2);
    my $n2 = $t2[0];
    my $s2 = $t2[1];
	
    if ($check) {
	my ($n11, $n12);
	my ($n21, $n22);
	my $dummy;
	($n11, $dummy, $n12) = ($n1 =~ /(.*)(\s+|\/)(1|2)/);
	($n21, $dummy, $n22) = ($n2 =~ /(.*)(\s+|\/)(1|2)/);
	
	die ">$n11< ne >$n21<" if ($n11 ne $n21);
	die ">$n12< == 1" unless ($n12 == 1);
	die ">$n22< == 2" unless ($n22 == 2);
    }


    if ($s1 eq 'Pass' && $s2 eq 'Pass') {
	&find_and_print($sh1, $oh1, $n1);
	&find_and_print($sh2, $oh2, $n2);
    } elsif ($s1 eq 'Pass') {
	&find_and_print($sh1, $ph1, $n1);
    } elsif ($s2 eq 'Pass') {
	&find_and_print($sh2, $ph2, $n2);
    }
}

close($fh1);
close($fh2);
close($sh1);
close($sh2);
close($oh1);
close($oh2);
close($ph1);
close($ph2);


print "Done!\n";


sub find_and_print {
    my ($ih, $oh, $name) = @_;

    while (defined(my $sl = <$ih>)) {
	chomp($sl);
	my @t = split(/\t/, $sl);
	my $n = $t[0];

	if (substr($n, 1) eq $name) {
	    print $oh "$sl\n";
	    $sl = <$ih>;
	    print $oh "$sl";
	    $sl = <$ih>;
	    print $oh "$sl";
	    $sl = <$ih>;
	    print $oh "$sl";
	    last;
	}
    }

}
