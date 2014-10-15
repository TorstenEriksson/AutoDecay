#!/usr/bin/perl -w
# ******************************************
# AutoDecay (version 5.06 for Perl)
# -- A utility to simplify the calculation of
# decay indices with the help of PAUP*
#
# Copyright (c) 2005 Torsten Eriksson
# torsten@utsteg.se
# This program is free software;
# you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# http://www.gnu.org/copyleft/gpl.html
# *******************************************
use strict;

	# -- For use with MacPerl
	# Uncomment the 7 lines below this one
	#if( $MacPerl::Version =~ /Application$/ ) {
	#        my( $cmdLine, @args );
	#        $cmdLine = &MacPerl::Ask( "Enter command line options:" );
	#        require "shellwords.pl";
	#        @args = &shellwords( $cmdLine );
	#        unshift( @ARGV, @args );
	#    }
if (not defined ($ARGV[0])) {
	die "Error: No option given. Use -s or -e or -h (help)\n";
	}
my($opt) = $ARGV[0]; # option
if (($opt ne "-e") and ($opt ne "-s")) {
	if ($opt ne "-h") {
		print "Error: Unknown option\n"; }
	print "Usage setup: AutoDecay -s tree-filename [paup command string within double quotes]\n";
	print " extraction: AutoDecay -e adc-filename tree_length\n";
	print "       help: Autodecay -h\n";
	exit(1);
	}
#--------------------------
if (not defined ($ARGV[1])) {
	die "Error: No file name given.\n";
	}
my($fil) = $ARGV[1]; # file name
my($version) = "5.06"; # version number

my($defaults);
if ($^O eq "linux" || $^O eq "darwin") {
	$defaults = "$ENV{'HOME'}/.AutoDecay";
	}
else {
	$defaults = "AD.def";
	}
my($adcsuffix) = "adc";
my($logsuffix) = "log";
my($tresuffix) = "tre";
my($lissuffix) = "dcy";
my $paup_default_cmds = "hsearch addseq=random nreps=100";
#
# Remove suffix, if present and fix filenames for output files
my($base);
$_ = $fil;
if (/\./) {
	/\.\w+$/;
	$base = $` . ".";
	}
	else {
	$base = $fil . "."; }
my($adcfile) = $base . $adcsuffix;
my($logfile) = $base . $logsuffix;
my($trefile) = $base . $tresuffix;
my($lisfile) = $base . $lissuffix;
my($cnt);
my($inputTree,$tree,$convTree);
my($nroftaxa,$ctcount,$ptr,$chr,$lastptr,$tp);
my(@ctree);

#-----------------------------------------------------------------------------
# ===========
#    SETUP
# ===========
if ($opt eq "-s") {
	my(@transNrs,@transTaxa);
	my($nr,$nam);
	my(@cdone);
	# ------------------------
	# Get PAUP* command string
	my($paupcmds) = $ARGV[2];
	if (not defined($paupcmds)) {
		# read defaults from file if string not supplied
		if (! -e $defaults) {
			warn "Warning: No paup* command found. Using defaults\n";
			$paupcmds = $paup_default_cmds;
			}
		else {
			open (DF,$defaults) || die "Error: Unable to open settings file\n";
			$paupcmds = <DF>;
			chomp ($paupcmds);
			close (DF) || warn "Warning: Could not close settings file\n";
			}
		}
	$_ = $paupcmds;
	/\w+/;
	my($searchtype) = $&;
	$paupcmds = $';
	# ----------------------------------
	# Get TREE from the input Nexus file
	my $block;
	my @treeblock_cmds;
	unless (get_block(\$block,"TREES",$fil)) { die "Error: No trees block found in file $fil"; }
	@treeblock_cmds = get_commands_array(\$block);
	my $translate = 0;
	my $transtable;
	my $tree_found = 0;
	foreach(@treeblock_cmds) {
		if (/^translate/i) {
			$translate = 1; # SET TO TRUE
			$transtable = $_;
			}
		elsif (/^\w*tree/i) {
			$tree_found = 1; # SET TO TRUE
			$inputTree = $_;
			}
		}
	if (! $tree_found) { die "Error: No tree found in file\n"; }
	$tree = nex_to_newick($inputTree); # remove nexus stuff from tree
	# ----------------------------------
	# Translate taxon names if a translation table present
	if ($translate) {
		my @tree_array;
		my @tt;
		# Put table in array
		while ($transtable =~ /(\d+)\s+([\w\d\.\-\+\_]+)/ig) {
			$tt[$1]=$2;
			}
		# put tree in array
		split_tree(\@tree_array,\$tree);
		# pass though tree array and exchange taxon numbers for names
		my $expecting = 0;
		for (my $i=0; $i<@tree_array;$i++) {
			$_ = $tree_array[$i];
			if ($expecting && /(\d+)/) {
				$tree_array[$i] = $tt[$1];
				$expecting = 0;
				}
			elsif (/[\(\,]/) { $expecting = 1; }
			else { $expecting = 0; }
			}
		# put tree back into a string (so that all the stuff below don't have to be rewritten)
		$tree = join("",@tree_array);
		}
	# Strip branch lengths from the tree
	$tree =~ s/:[\d.]+//g && warn "Information: Branch lengths stripped from the input tree\n";
	#---------------------------------------
	# Create CONSTRAINT trees for each clade
	$_ = $tree;
	$nroftaxa = tr/\,//; # count number of taxa in tree by counting commas + 1
	$nroftaxa++;
	$ctcount = 0;
	$ptr = 0;
	while ($ptr <= length($tree)) {
	   $chr = substr($tree,$ptr,1); # get a character
	   if ($chr eq "(") {
	      # start a new tree
	      $ctree[$ctcount] = "((";
	      $cdone[$ctcount] = 0;
	      ++$ctcount; }
	   elsif ($chr eq ")") {
	      # end a tree
 	      $lastptr = $ctcount;
	      --$lastptr; #point to last tree
	      while ($cdone[$lastptr] == 1) { # step past possible done trees
	         --$lastptr; }
	      $ctree[$lastptr] .= "))";
	      $cdone[$lastptr] = 1;
	      }
   	   else {
          # copy character to all active trees
          $tp = 0;
          while ($tp < $ctcount) {
    	     if ($cdone[$tp] == 0) {
                $ctree[$tp] .= $chr;
        	    }
             ++$tp;
          }
       }
       ++$ptr;
    }
    --$ctcount;
    # -------------------------------
	# Scan and mark meaningless trees
	# (those with less than two taxa inside or outside)
	$cnt = 0;
	while ($cnt < @ctree) {
		$_ = $ctree[$cnt];
		$tp = tr/\,//;
		if (($tp < 1) or ($tp >= $nroftaxa-2)) {
			$cdone[$cnt] = 2; # mark as not usable
		}
		$cnt++;
	}
	# WRITE paup data to adc file
	open(AF,">$adcfile") || die "Error: Could not open new adc file\n";
	print AF "#NEXUS\n";
	print AF "BEGIN PAUP;\n";
	print AF "log file='$logfile' start replace;\n";
	print AF "set WarnRedef=no;\n";
	chomp;
	print AF "[!AutoDecay ($version) command file]\n";
	print AF "[!Input tree: $inputTree]\n";
	print AF "[!From file: $fil]\n";
	if ($translate) {
		print AF "[!Translated input tree: $tree]\n\n";
		}
	print AF "set autoclose;\n";
	$cnt = 0;
	$tp = 1;
	while ($cnt <= $ctcount) {
		if (not $cdone[$cnt] == 2) {
			print AF "[!---> Constraint tree # $tp]\n";
			print AF "Constraints AD$tp = $ctree[$cnt];\n";
			print AF "$searchtype constraints=AD$tp enforce converse $paupcmds;\n";
			$tp++;
		}
		$cnt++;
	}
	print AF "log stop;\n";
	print AF "END;\n";

	# write AutoDecay data to adc file
	print AF "\nBegin AutoDecay;\n";
	$cnt = 0;
	while ($cnt <= $ctcount) {
		if ($cdone[$cnt] == 2) {
			print AF "$cnt "; # save info of what trees were not used
		}
		$cnt++;
	}
	print AF "\n";
	print AF "$tree\n"; #save (possibly converted) tree
	$cnt = 0;
	while ($cnt <= $ctcount) {
		print AF "$ctree[$cnt]\n";
		$cnt++;
	}
	print AF "end;\n";
	close(AF) || warn "Warning: Could not close adc file\n";

	#
	exit(0);
	}
#-----------------------------------------------------------------------------
# ================
#    EXTRACTION
# ================
else ## ($opt eq "-e")
	{
	if (not defined ($ARGV[2])) {
		die "Error: Please input optimal tree length\n";
		}
	my($optimal) = $ARGV[2]; # score of optimal tree
	my(@cscore, @dvalue);
	my(@delTrees);
	my($ix);
	my(@ignore);
	my(@nodePtr);
	my($newTree);
	# ----------------------------------
	# Get constraint trees and other data from adc file
	# delTrees inputTree trees
	my $block;
	my @constraint_trees;
	unless (get_block(\$block,"AUTODECAY",$adcfile)) { die "Error: No AutoDecay block found in file $adcfile"; }
	@ctree = split(/\s+/,$block);
	shift(@ctree); # remove initial begin autodecay block command
	shift(@ctree); # remove initial begin autodecay block command
	pop(@ctree);   # remove end command
	$cnt=0;
	foreach(@ctree) {
		if (/^\(/) { last; }
		$cnt++;
		}
	@delTrees = splice(@ctree,0,$cnt); # LIst of trees to be ignored
	$inputTree = shift(@ctree); # The first tree is original input tree
	# -------------------------------
	# Adjust info about ignored trees
	for (my $i=0;$i<@ctree;$i++) {
		$ignore[$i]=0; # set all to be used by default
		}
	foreach(@delTrees) {
		$ignore[$_] = 1; # set tree to be ignored
		}
	# ----------------------------
	# Read tree info from log file
	#
	open (LF,$logfile) || die "Error: Could not open log file\n";
	$cnt=0;
	$ix=0;
	while(<LF>) {
		chomp;
		if ((/Score of best tree/) or (/Length of shortest tree/)) {
			# A completed analysis found
			/[\d\.]+$/;
			$cscore[$cnt] = $&; # tree score
			while ($ignore[$ix]) {$ix++;} # step past those to be ignored
			$dvalue[$ix] = $cscore[$cnt] - $optimal; # compute decay value
			$cnt++;
			$ix++;
			}
		}
	close (LF) || warn "Warning: Could not close log file\n";
	# Check that anything was found
	if (not defined($cscore[0])) {
		die "Error: No analysis results found in log file\n";
		}
	if (@cscore != @ctree - @delTrees) {
		die "Error: Tree number mismatch in log file\nlength(@cscore) length(@ctree) length(@delTrees)\n"
		}
	# ---------------------------------
	# Write decay value list (adjust for deleted trees delTrees)
	#
	open (LF,">$lisfile") || die "Error: Could not open output file for writing.\n";
	$cnt = 0;
	while ($cnt < @ctree) {
		if (not ($ignore[$cnt])) {
			print LF "D=$dvalue[$cnt] Node=$ctree[$cnt]\n";
			}
		$cnt++;
		}
	close (LF) || die "Error: Could not close output file 1\n";
	# ----------------------------------------------
	# Write tree with values as internal node labels
	open (LF,">$trefile") || die "Error: Could not open output file for writing.\n";
	print LF "#NEXUS\n";
	$_ = `date`;
	chomp;
	print LF "[! Tree file created by AutoDecay ($version) at $_]\n";
	print LF "[! Decay indices are indicated as internal node labels in the tree]\n";
	print LF "Begin trees;\n";
	print LF "\tTree Decay = ";
	$ix = 0;
	$ptr = 0;
	$newTree = "";
	while ($ptr < length($inputTree)) {
		$_ = substr($inputTree,$ptr,1); # get a character
		if (/\(/) {
			$newTree .= "(";
			push(@nodePtr,$ix);
			$ix++;
			}
		elsif (/\)/) {
			$tp = pop(@nodePtr);
			if (not $ignore[$tp]) {
				$newTree .= ")'d " . $dvalue[$tp] . "'";
				}
			else {$newTree .= ")";}
			}
		else {
			$newTree .= $_;
			}
		$ptr++;
		}
	print LF "$newTree;\n";
	print LF "End;\n";
	close (LF) || die "Error: Could not close output file 2\n";

	exit(0);
	}

# ==============================================================
#   ---- SUBROUTINES THAT DEAL WITH NEXUS FILE AND TREES -----
# ==============================================================

# ====================================
# Get a specific BLOCK from nexus data
# ====================================
sub get_block {
	# Input:	Ref to string where the block is to be stored
	#			Name of block
	#			Name of nexus file
	# Output:	Block (including initial "Begin xxx;" and ending "end;") is
	#			   stored in the referenced string.
	# Return:	1 if OK, 0 if not found or other error
	my($theblock,$block_name,$nexus_file) = @_;
	my(@blocks);
	get_blocks_array(\@blocks,$nexus_file) || return 0;
	foreach(@blocks) {
		if (/^Begin\s+$block_name/i) {
			$$theblock = $_;
			return 1;
			}
		}
	return 0; # No such block found
	}
# ================================
# Get all BLOCKs from a nexus file
# ================================
sub get_blocks_array {
	# Input:	Ref to an array where blocks are to be stored
	#			Nexus file name
	# Output:	Blocks are stored in the referenced array:
	#			-- One block in each element including the initial
	#			   "begin xxx;" and ending "end;" commands.
	# Return:	1 if OK, 0 if an error occurred
	#
	my($blocks,$nexfile) = @_;
	my($nexdata);
	get_clean_nexus(\$nexdata,$nexfile) || return 0;
	@$blocks = split(/end\s*;|endblock\s*;/i,$nexdata);
	$$blocks[0] =~ s/^\s*#NEXUS//i;  # remove starting statement if any
	# Reset the end command to each block
	my($nr_of_blocks) = scalar @$blocks;
	my($i);
	for ($i=0; $i < $nr_of_blocks ; $i++) {
		$$blocks[$i] =~ s/^\s*//;	# remove leading whitespace
		if (! $$blocks[$i] eq "") {
			$$blocks[$i] .= " End;";		# reset end command
			}
		}
	# remove any empty list items at end
	if ($$blocks[$nr_of_blocks-1] eq "") {
		pop(@$blocks);
		}
	return 1;
	}
# ==================================
# Get cleaned data from a nexus file
# ==================================
sub get_clean_nexus {
	# Input: 	Ref to string in which to store data
	#			file name
	# Output:	Cleaned data (without nexus comments) stored in referenced string
	# Return:	1 when OK, or 0 on error
	# Data treatment:
	#   * Nexus style comments are removed (comments in comments are not handled)
	#   * Initial and trailing whitespace is removed
	# --------------
	# Get nexus data
	# --------------
	my($data,$nexfile) = @_;
	get_nexus($data,$nexfile) || return 0;
	# -----------------------------
	# Strip comments and whitespace
	# -----------------------------
	$$data =~ s/\[.*?\]//g;	# remove nexus comments
	$$data =~ s/^\s*//;		# remove leading whitespace
	$$data =~ s/\s*$//;		# remove trailing whitespace
	return 1;
	}
# ==========================
# Get data from a nexus file
# ==========================
sub get_nexus {
	# Input:	Ref to string in which to store data
	# 			file name
	# Output:	Data stored without linebreaks in referenced string
	# Return:	returns 1 when OK, or 0 on error
	# Data treatment:
	#   * Line breaks (newline and/or carriage-returns) are substituted into spaces
	# ------------------
	# Open and read file
	# ------------------
	my($nexdata,$nexfile) = @_;
	open (NF,$nexfile) || return 0;
	my(@nex) = <NF>;
	close (NF);
	$$nexdata = join(" ",@nex);
	# -----------------
	# Strip line breaks
	# -----------------
	$$nexdata =~ s/[\n\r]/ /g;	# substitute line breaks into spaces
	return 1;
	}
# ==========================================
# Get commands from a nexus block into array
# ==========================================
sub get_commands_array {
	# Input:	Ref to string with a block
	# Output:	Array of commands
	# Note that the ending ";" is removed in all commands
	# Data treatment: leading and trailing whitespace is removed for each command
	my($block) = @_;
	my(@cmds) = split (/;/,$$block);
	foreach (@cmds) {
		s/^\s*//; # remove leading whitespace
		s/\s*$//; # remove trailing whitespace
		}
	my($nr_of_cmds) = scalar @cmds;
	# remove any empty list items at end
	if ($cmds[$nr_of_cmds-1] eq "") {
		pop(@cmds);
		}
	return @cmds;
	}
# =====================================
# Transform Nexus tree to Newick format
# =====================================
# Removes nexus stuff from a proper tree
# This does not strictly follow a Newick tree since it lacks a ";" at end
sub nex_to_newick {
	# Input:	a nexus tree command
	# Output:	cleansed tree
	$_ = $_[0];
	s/\[.*?\]//g;	# remove any embedded nexus comments
	/\(.*\)/; # match tree
	return $&;
	}
# =================================
# Split a newick tree into an array
# =================================
sub split_tree {
	# Input:	Ref to array where tree is to be stored
	#			Ref to tree in Newick format
	# Output:	Tree stored in the referenced array
	# Return:	1 if OK, 0 if an error occurred
	my($split_tree,$newick_tree) = @_;
	my($pointer)=0;
	my($token);
	while() {
		$token = nxt_token($newick_tree,\$pointer);
		if ($token eq "") { last ; }
		push(@$split_tree,$token);
		}
	return 1;
	}

sub nxt_token {
	# Input:	Ref to a string
	#			Ref to a pointer of where to commence
	# Output:	Pointer updated
	# Return:	A token or undef if at end of string
	my($data,$ptr) = @_;
	my($tmp);
	my($token)="";
	if ($$ptr < length($$data)) {
		while ($$ptr < length($$data)) {
			$tmp = substr($$data,$$ptr,1); # get a character
			if ($tmp =~ /[\(\),:;]/) {
				if ($token eq "") {
					$$ptr++;
					return $tmp;
					}
				else {
					return $token;
					}
				}
			elsif (($token ne "") and ($tmp =~ /\s/)) {
				$$ptr++;
				return $token;
				}
			elsif (not $tmp =~ /\s/) {
				$token .= $tmp;
				$$ptr++;
				}
			else {
				$$ptr++;
				}
			}
		if (not $token eq "") { return $token }
		}
	}
