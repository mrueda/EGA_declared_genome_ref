#!/usr/bin/env perl
#
#   Script to parse EGAZ-XML files
#
#   Last Modified; Jun/30/2021
#
#   Version 0.0.1
#
#   Copyright (C) 2021 Manuel Rueda (manuel.rueda@crg.eu)
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, see <https://www.gnu.org/licenses/>.
#
#   If this program helps you in your research, please cite

use strict;
use warnings;
use autodie;
use feature qw(say);
use Data::Dumper;

# ==========
# BACKGROUND
# ==========
# EGA submitters are not enforced to use a given nomenclature to declare reference genomes.
# The declared reference genome information is inside EGAZ XML files.
#
# ========
# STRATEGY
# ========
# The XMLs can contain the actual genome name, synonyms, flavours, or whatever.
# Such info can be written in > 1 XML element so it's mandatory to parse the XML.
#
# Before parsing anything, we carried out an exhaustive exploration of the XML contents and found:
#
# 1 - The genome reference could be declared in these two elements:
#  <ANALYSIS_TYPE>
#      <SEQUENCE_VARIATION>
#        <ASSEMBLY>
#          <STANDARD refname="GRCh37"/>   <==================================================== HERE
#        </ASSEMBLY>
#        <SEQUENCE accession="accession="GL000207.1" label="GL000207.1"/> <==================== HERE
#
# 2 - A significant number of XML 56% (63927 out 113852) belong to the project DDD (https://ega-archive.org/studies/EGAS00001000775)
#     From those, 57131 contain  the pattern 'URL/'
#     <SEQUENCE accession="" label="2,URL=/lustre/scratch113/projects/ddd/ref_genome/hs37d5/2,assembly=hs37d5,length=243199373"/>
#
# 3 - According to this, we created an algorithm where we fetch the declared genome sequentially,
#     depending on the presence/absance of fields:
#
#     * First use <STANDARD refname> to fetch the declared genome
#     * If not found then (INFT) we parse 'ref_genome' 
#     * INFT we check <SEQUENCE label>
#     * INFT we check <SEQUENCE accession>
#
#     The above search is perfomed line by line, thus, the first successful match "wins" ('label' can win over 'ref_genome').
#     Note that we use a dictionary to minimize the number of possible answers.
#
# ==============
# IMPLEMENTATION
# ==============
#
# There are many ways to get attribute values from XML (e.g., xmllint, etc):
# https://stackoverflow.com/questions/15461737/how-to-execute-xpath-one-liners-from-shell
# As our task is simple, here we have chosen Perl (usually installed by default in all Linux distros).
# Note that we'll parse the XML elements directly (w/o using XML parsers).

my $DEFAULT = 'NA';    # global variable
my %genome  = (        # global hash
    GRCh37 => 'GRCh37',
    GRCh38 => 'GRCh38',
    hg19   => 'GRCh37',
    hs37d5 => 'hs37d5',
    hg17   => 'hg17'
);

# We declare flavours of GRCh37
my @thirty_sevens = qw (
  GRCh37.p13
  GRCh37.p1
  GRCh37.decoy
  GCA_000001405.1
  GCA_000001405.5
  GL000207.1
  CM000683.11
  NC_000001.10
);

# We declare flavours of GRCh38 (note that a regex will also work - see below) 
my @thirty_eights = qw (
  CM000663.2
  CM000682.2
  CM000666.2
  CM000684.2
);

# Adding the flavours to %genome
@genome{@thirty_sevens} = ('GRCh37') x scalar @thirty_sevens;
@genome{@thirty_eights} = ('GRCh38') x scalar @thirty_eights;
#print Dumper \%genome;

# Defining a few more variables
my $ref_genome = $DEFAULT;
my $last       = 0;
my $pattern1   = 'STANDARD ';    # trailing space matters
my $pattern2   = 'SEQUENCE ';    # trailing space matters
my @patterns = ( 'ref_genome', 'label', 'accession' );    # order matters

#####################
# XML PARSING START #
#####################

XML: while (<>) {
    next unless /$pattern1|$pattern2/;
    my %data = ();
    if (/$pattern1/) {
        next if /HOMO/i;    # Some folks decided to write 'HOMO SAPIENS'
        ( $ref_genome, $last ) = assign_ref_gen( 'refname', $_, 0 );
        $data{$pattern1} = $ref_genome;
        last if $last;
    }
    elsif (/$pattern2/) {
        for my $pattern (@patterns) {
            my $switch = $pattern eq 'ref_genome' ? 1 : 0;
            ( $ref_genome, $last ) = assign_ref_gen( $pattern, $_, $switch );
            $data{$pattern} = $ref_genome;
            last XML if $last;
        }
    }

    #print Dumper \%data;
}
say $ref_genome;

###################
# XML PARSING END #
###################

sub assign_ref_gen {

    my ( $pattern, $line, $switch ) = @_;
    my $out   = '';
    my $last  = 0;
    my $regex = $switch ? qq{\/$pattern\/(\\w+)\/} : qq{$pattern="(\\S+)"};
    $line =~ /$regex/;
    my $ref = $1 // $DEFAULT;
    $ref = ( $ref =~ /CM000\d{3}.1/ ) ? 'GRCh37' : $ref;
    if ( exists $genome{$ref} ) {
        $out  = $genome{$ref};
        $last = 1;
    }
    return ( $out, $last );
}
