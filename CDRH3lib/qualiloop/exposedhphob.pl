#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    exposedhphob
#   File:       exposedhphob.pl
#   
#   Version:    V1.0
#   Date:       20.11.20
#   Function:   Calculate a quality score for exposure of hydrophobic/
#               hydrophilic residues
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2020
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
use strict;

# Check environment variables
if(!defined($ENV{'DATADIR'}))
{
    print STDERR "You must set DATADIR to the BiopLib/BiopTools data directory";
    exit 1;
}

# Set the hydrophobicity file
my $HPhobFile = $ENV{'DATADIR'} . "/consensus.hpb";
$HPhobFile = $::hphob if(defined($::hphob));

# Check command line
UsageDie($HPhobFile) if((scalar(@ARGV) != 3) || defined($::h));

# Take parameters from command line
my $startRes = shift(@ARGV);
my $stopRes  = shift(@ARGV);
my $pdbFile  = shift(@ARGV);

# Read the hydrophobicity file
my ($minHPhob, $maxHPhob, %hphobs) = ReadHPhobs($HPhobFile);

# Normalize such that the scores go from -1 to +1.
%hphobs = NormalizeHPhobs(%hphobs);

# Test code to replace RunPdbsolv
#my ($aResIds, $aSequence, $aAccess) = CreateTest();

# Calculate accessibility and extract desired region
my ($aResIds, $aSequence, $aAccess) = RunPdbsolv($pdbFile, $startRes, $stopRes);

# Calculate quality scores and print them
CalculateQualityAndPrint($aResIds, $aSequence, $aAccess, %hphobs);

#*************************************************************************
#> ($aResIds, $aSequence, $aAccess) = CreateTest();
#  ------------------------------------------------
#  Returns: $aResIds    reference to array of residue IDs
#           $aSequence  reference to array of amino acid names
#           $aAccess    reference to array of acccessibilities
#
#  Test routine that creates a simple test with buried and exposed
#  residues that are very hydrophobic, slightly hydrophobic,
#  slightly hydrophilic, very hydrophilic.
#
#  20.11.20 Original   By: ACRM
sub CreateTest
{
    my @resIds   = qw/H1 H2 H3 H4 H5 H6 H7 H7/;
    my @sequence = qw/ILE ILE CYS CYS PRO PRO ARG ARG/;
    my @access   = qw/0.01 1.0 0.01 1.0 0.01 1.0 0.01 1.0/;
        
    return(\@resIds, \@sequence, \@access);
}

#*************************************************************************
#> void CalculateQualityAndPrint($aResIds, $aSequence, $aAccess, %hphobs);
#  -----------------------------------------------------------------------
#  Input:   $aResIds    reference to array of residue IDs
#           $aSequence  reference to array of amino acid names
#           $aAccess    reference to array of acccessibilities
#           %hphobs     hash of hydrophobicities indexed by AA name
#
#  Calculate the quality scores and print them
#
#      /  1 + (H * (1-A))   if (H<0)
#  Q = |
#      \  1 - (H * A)       otherwise
#
#  where A = accessibility (0..1) and H = hydrophobicity (-1..+1)
#
#  This has the desired property of being close to 1 if we have a buried
#  hydrophobic, an accessible hydrophilic, or if the hydrophobicity is
#  close to 0.
#
#  i.e. if Q=f(A,H)
#
#  Hydrophilics    Hydrophobics
#  f(0,  0) = 1    f(0, 0) = 1
#  f(0, -1) = 0    f(0, 1) = 1
#  f(1,  0) = 1    f(1, 0) = 1
#  f(1, -1) = 1    f(1, 1) = 0
#
#  20.11.20 Original   By: ACRM
sub CalculateQualityAndPrint
{
    my($aResIds, $aSequence, $aAccess, %hphobs) = @_;

    my $nRes        = scalar(@$aResIds);
    my $totalQScore = 0.0;
    
    for(my $i=0; $i<$nRes; $i++)
    {
        my $qualityScore;
        my $resHPhob = $hphobs{$$aSequence[$i]};

        if($resHPhob < 0)
        {
            $qualityScore = 1 + ($resHPhob * (1 - $$aAccess[$i]));
        }
        else
        {
            $qualityScore = 1 - ($resHPhob * $$aAccess[$i]);
        }
        
        printf "%-6s %3s %5.3f\n",
            $$aResIds[$i], $$aSequence[$i], $qualityScore;

        $totalQScore += $qualityScore;
    }
    printf "Total: %.3f\n", $totalQScore;
    printf "Mean:  %.3f\n", $totalQScore / $nRes;
}


#*************************************************************************
#> ($aResIds, $aSequence, $aAccess) =
#        RunPdbsolv($pdbFile, $startRes, $stopRes);
#  ------------------------------------------------
#  Input:   $pdbFile    PDB filename
#           $startRes   Res id of first residue of interest
#           $stopRes    Res id of last residue of interest
#  Returns: $aResIds    reference to array of residue IDs
#           $aSequence  reference to array of amino acid names
#           $aAccess    reference to array of acccessibilities
#
#  Test routine that creates a simple test with buried and exposed
#  residues that are very hydrophobic, slightly hydrophobic,
#  slightly hydrophilic, very hydrophilic.
#
#  20.11.20 Original   By: ACRM
sub RunPdbsolv
{
    my($pdbFile, $startRes, $stopRes) = @_;
    my @sequence = ();
    my @access   = ();
    my @resIds   = ();
    
    my $tfile = "/var/tmp/exhpb_" . $$ . time();

    `pdbsolv -r $tfile $pdbFile`;
    if(open(my $fp, '<', $tfile))
    {
        my $inRange = 0;
        while(<$fp>)
        {
            my $resId     = substr($_, 8, 7);
            $resId        =~ s/\s//g;
            my $resAccess = substr($_, 30, 7);
            my $aa        = substr($_, 17, 3);
            $inRange = 1 if($resId eq $startRes);
            if($inRange)
            {
                push @resIds,   $resId;
                push @sequence, $aa;
                push @access,   $resAccess / 100.0;
            }
            $inRange = 0 if($resId eq $stopRes);
        }
        close $fp;
    }
    else
    {
        print STDERR "Error: Can't open solvent accessibility file ($tfile)\n";
        exit 1;
    }
    unlink $tfile;
    
    return(\@resIds, \@sequence, \@access);
}


#*************************************************************************
#> %hphobs = NormalizeHPhobs(%hphobs);
#  -----------------------------------
#  Input:   %hphobs     hash of hydrophobities indexed by amino acid name
#  Returns: %hphobs     hash of hydrophobities indexed by amino acid name
#
#  Normalizes the hydrophobities to a range of -1 to +1
#
#  20.11.20 Original   By: ACRM
sub NormalizeHPhobs
{
    my(%hphobs) = @_;
    return(SignNormalizeHPhobs(%hphobs));
#    return(FullNormalizeHPhobs(%hphobs));
}

#*************************************************************************
#> %hphobs = SignNormalizeHPhobs(%hphobs);
#  ---------------------------------------
#  Input:   %hphobs     hash of hydrophobities indexed by amino acid name
#  Returns: %hphobs     hash of hydrophobities indexed by amino acid name
#
#  Normalizes the hydrophobities to a range of -1 to +1 Normalizes the
#  negative range (-1..0) and the positive range (0..+1) separately
#
#  20.11.20 Original   By: ACRM
sub SignNormalizeHPhobs
{
    my(%hphobs) = @_;
    
    my $HPhobSpan = $maxHPhob - $minHPhob;

    foreach my $aa (sort keys %hphobs)
    {
        my $value = $hphobs{$aa};
        if($value < 0.0)
        {
            $value /= $minHPhob;
            $value *= -1;
        }
        else
        {
            $value /= $maxHPhob;
        }
        $hphobs{$aa} = $value;
    }
    return(%hphobs);
}

#*************************************************************************
#> %hphobs = FullNormalizeHPhobs(%hphobs);
#  -----------------------------------
#  Input:   %hphobs     hash of hydrophobities indexed by amino acid name
#  Returns: %hphobs     hash of hydrophobities indexed by amino acid name
#
#  Normalizes the hydrophobities to a range of -1 to +1 over the whole
#  input range
#
#  20.11.20 Original   By: ACRM
sub FullNormalizeHPhobs
{
    my(%hphobs) = @_;
    
    my $HPhobSpan = $maxHPhob - $minHPhob;

    foreach my $aa (sort keys %hphobs)
    {
        my $value = $hphobs{$aa};
        $value -= $minHPhob;
        $value /= $HPhobSpan;
        $value *= 2.0;
        $value -= 1.0;
        $hphobs{$aa} = $value;
    }
    return(%hphobs);
}


#*************************************************************************
#> ($minHPhob, $maxHPhob, %hphobs) = ReadHPhobs($HPhobFile);
#  ---------------------------------------------------------
#  Input:   $HPhobFile  Hydrophobicity filename
#  Returns: $minHPhob   minimum hydrophobicity value
#           $maxHPhob   maximum hydrophobicity value
#           %hphobs     hash of hydrophobities indexed by amino acid name
#
#  Reads the specified hydrophobicity file
#
#  20.11.20 Original   By: ACRM
sub ReadHPhobs
{
    my($filename) = @_;

    my %hphobs = ();
    my $minHPhob =  1000.0;
    my $maxHPhob = -1000.0;
    if(open(my $fp, '<', $filename))
    {
        <$fp>;  # Skip title line
        while(<$fp>)
        {
            chomp;
            my @fields = split(/\s+/, $_);
            $hphobs{$fields[0]} = $fields[1];
            $minHPhob = $fields[1] if($fields[1] < $minHPhob);
            $maxHPhob = $fields[1] if($fields[1] > $maxHPhob);
        }
        close($fp);
    }
    else
    {
        print STDERR "Error: Can't read hydrophobicity file - $filename\n";
        exit 1;
    }
    return($minHPhob, $maxHPhob, %hphobs);
}

#*************************************************************************
#> UsageDie()
#  ----------
#  Prints a usage message and exits
#
#  20.11.20 Original   By: ACRM
sub UsageDie
{
    my($HPhobFile) = @_;
    
    print <<__EOF;

exposedhphob.pl V1.0 (c) 2020 UCL, Prof. Andrew C.R. Martin.

Usage: exposedhbob.pl [-hphob=hphobfile] startres stopres file.pdb
       -hphob   - specify hydrophobicity file [$HPhobFile]
       startres - first residue ID in the format [c[.]]nnn[i]
       stopres  - last residue ID in the format [c[.]]nnn[i]
       file.pdb - PDB file
    
Outputs a 'quality' score (0-1) for each residue in the specified range.
This is calculated as

      /  1 + (H * (1-A))   if (H<0)
  Q = |
      \\  1 - (H * A)       otherwise

where A = accessibility (0..1) and H = hdrophobicity (-1..+1)

The total and average scores are also output.

The calculation favours buried hydrophobics and exposed hydrophilics;
it penalizes exposed hydrophobics and buried hydrophilics.

__EOF
    exit 0;
}

