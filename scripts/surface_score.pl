#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables

 sub cal_fractional_exposed($);        # get the fractional exposed area, check the paper
 sub convert_dssp_to_three_ss($); # convert a secondary structure of dssp into helix(G,H,I), strand(E,B), and loop(others).

  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args

	  print"This script process the secondary structure parsed by DSSP. Use the fraction area of exposed nonpolar residues as the quality score!\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";
	  print "0. Dir of secondary structure parsed by dssp.!\n";
	  print "1. Dir of output\n";


	  print "\nFor example:\n";

           print "\n**************** ab initio for validation *****************\n";
          print "perl $0 ../dssp_parsed_DB ../1_calculated_scores/feature_0_surface_DB\n";
	  exit(0);
	}
############### the copy right for each script #############################

############################################################################

#################################################################################################################################
my(@aanames) = ('A', 'B', 'C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X', 'Y', 'Z');
my(@accth) = (106, 160, 135, 163, 194, 197, 84, 184, 169, 205, 164, 188, 157, 136, 198, 248, 130, 142, 142, 227, 180, 222, 196);
my(%accth2)=();
my($i);
for ($i = 0; $i<=$#aanames; $i++) {
	$accth2 {$aanames[$i]} = $accth[$i];
}
##################################################################################################################################
my(%aa_polar)=();              # tell whether a aminoacid is polar or not. If is polar, value is 1, else 0
$aa_polar{'A'} = 0;
#$aa_polar{'B'} = 0;
$aa_polar{'C'} = 0;
$aa_polar{'D'} = 1;
$aa_polar{'E'} = 1;
$aa_polar{'F'} = 0;
$aa_polar{'G'} = 0;
$aa_polar{'H'} = 1;
$aa_polar{'I'} = 0;
$aa_polar{'K'} = 1;
$aa_polar{'L'} = 0;
$aa_polar{'M'} = 0;
$aa_polar{'N'} = 1;
$aa_polar{'P'} = 0;
$aa_polar{'Q'} = 1;
$aa_polar{'R'} = 1;
$aa_polar{'S'} = 1;
$aa_polar{'T'} = 1;
$aa_polar{'V'} = 0;
$aa_polar{'W'} = 0;
#$aa_polar{'X'} = 0;
$aa_polar{'Y'} = 1;
#$aa_polar{'Z'} = 0;
##############standard Amino Acids (3 letter <-> 1 letter)#######
my(%amino)=();
$amino{"ALA"} = 'A';
$amino{"CYS"} = 'C';
$amino{"ASP"} = 'D';
$amino{"GLU"} = 'E';
$amino{"PHE"} = 'F';
$amino{"GLY"} = 'G';
$amino{"HIS"} = 'H';
$amino{"ILE"} = 'I';
$amino{"LYS"} = 'K';
$amino{"LEU"} = 'L';
$amino{"MET"} = 'M';
$amino{"ASN"} = 'N';
$amino{"PRO"} = 'P';
$amino{"GLN"} = 'Q';
$amino{"ARG"} = 'R';
$amino{"SER"} = 'S';
$amino{"THR"} = 'T';
$amino{"VAL"} = 'V';
$amino{"TRP"} = 'W';
$amino{"TYR"} = 'Y';
###################################################################



 my($addr_input)=$ARGV[0];
 my($addr_output)=$ARGV[1];

 -s $addr_input || die "cannot open input $addr_input\n";





 my(%hash)=();
 my($path_target,$search_name,$path_model1,$path_model2,$path_search,$search,$path_out,$name,$path_write,$file,$return_val,$seq1,$seq2,$quality);
 my(@tem_split,@targets,@files,@searches);
 my($IN,$line,$OUT);
###########################################################################################################

 my($target);
 ############# read the dssp parsed secondary structure for each model ###############
	$path_target=$addr_input;          # the target folder

    $path_out = $addr_output;      # the output file
	$OUT=new FileHandle ">$path_out";
	defined($OUT) || die "Cannot open output $path_out\n";
    opendir(DIR,"$path_target");
	@targets=readdir(DIR);
	foreach $target (@targets)
	{
		 if($target eq '.' || $target eq '..')
		 {
			 next;
		 }
		 $path_model1=$path_target."/".$target;             # this is the dssp parsed secondary structure.
         @tem_split=split(/\./,$target);         # try to get the target name, remove .dssp_parsed
		 $name=$tem_split[0];
		 for($i=1;$i<@tem_split-1;$i++)
		 {
			 $name.=".".$tem_split[$i];
		 }

		 ######### here can be changed! For quality calculation     ###########
		 $quality = cal_fractional_exposed($path_model1);
		 if($quality == -1)
		 {# this model is not correctly parsed by dssp or something is wrong, we skip to give score to these models.
			 next;
		 }

                 $quality = 1-$quality;

		 print $OUT $name."\t".$quality."\n";
	}
	$OUT->close();

 sub convert_dssp_to_three_ss($)
 {
	 my($cha)=@_;
	 if($cha eq "G" || $cha eq "H" || $cha eq "I")
	 {# helix
		 return "H";
	 }
	 elsif($cha eq "E" || $cha eq "B")
	 {#strand
		 return "E";
	 }
	 else
	 {#coil
		 return "C";
	 }
 }

 sub cal_fractional_exposed($)
 {# read dssp processed result and parse the exposed area for each residue
	 my($input)=@_;
	 my($IN,$line);
     my(@tem_split);
	 my(@aa)=();                        # the amino acid
	 my(@exposed)=();                   # exposed area for each residue

	 $IN=new FileHandle "$input";
	 if(defined($line=<$IN>))
	 {
		 # this is for chain ID
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is for total number
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is for amino acid
		 chomp($line);
		 @aa=split(/\s+/,$line);       # this is the amino acid
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is for index
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is secondary structure

	 }
	 if(defined($line=<$IN>))
	 {
		 # this is BP1

	 }
	 if(defined($line=<$IN>))
	 {
		 # this is BP2

	 }
	 if(defined($line=<$IN>))
	 {
		 # this is solvent
		 chomp($line);
		 @exposed=split(/\s+/,$line);  # this is the exposed area
	 }
	 $IN->close();
     my($ex_non)=0;           # the total exposed area of nonpolar residue
	 my($ex_all)=0;           # the total exposed area of all residues
	 my($index)=scalar(@aa);
	 my($i);
	 for($i=0;$i<$index;$i++)
	 {
                 if(not exists $aa_polar{$aa[$i]}) {next;}
		 $ex_all+=$exposed[$i];
		 if($aa_polar{$aa[$i]} == 0)
		 {# this is a nonpolar residue
			 $ex_non+=$exposed[$i];
		 }
	 }
	 if($ex_all == 0)
	 {
######### I only get few models which dssp cannot parse, because the pdb format is not correct! To make it simple, I give these model -1, it may influence a little bit about the final performance #####
		 print "Warning! Error, check $input, no residue inside!\n";
		 return -1;
		 #exit(0);
	 }
	 my($quality)=$ex_non/$ex_all;
	 return $quality;

 }
