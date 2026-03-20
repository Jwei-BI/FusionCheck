#!/usr/bin/perl
use strict;
use Cwd;
use File::Spec;
use Getopt::Long;
use feature qw(say);
use Time::HiRes qw(gettimeofday tv_interval);

# Default parameter values
my $blat        = "/software/blat";
my $blat_args   = "-stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -t=dna -q=dna -out=psl";
my $hg19        = "/ref/hg19.2bit";
my $samtools    = "/software/samtools/samtools-1.9/samtools";
my $tabix       = "/software/samtools/samtools-1.9/htslib-1.9/tabix";
my $refgene     = "/ref/refGene.txt.sort.gz";
my $extend      = 10;
my $match_min   = 30;
my $out_dir     = getcwd();
# my $bam = "/sds1/Project/commercial_pro/Project_ctdna/ctDNA_1123_SJ/20240513/BJ24CM001538T7/analysis/BJ24CM001538T7/BJ24CM001538T7_T/Map/BJ24CM001538T7_T.clipped.bam";
# my @breakpoint = qw(chr20:43961871 chr6:117650683);
# my $prefix = "BJ24CM001538T7";
my $regions;
my $bam_file;
my $prefix = "fusionCheck";
my $help;
# Parse command line parameters
GetOptions(
    'bam=s'         => \$bam_file,
    'regions=s'     => \$regions,
    'extension=i'   => \$extend,
    'min_match=i'   => \$match_min,
    'blat=s'        => \$blat,
    'samtools=s'    => \$samtools,
    'tabix=s'       => \$tabix,
    'hg19=s'        => \$hg19,
    'refgene=s'     => \$refgene,
    'prefix=s'      => \$prefix,
    'out_dir=s'     => \$out_dir,
    'help'          => \$help,
) or die "Incorrect usage!\n";

# Display help information
if ($help || !$bam_file || !$regions) {
    print <<HELP;
Usage: $0 --bam <bam_file> --regions <"region1 region2"> [--extension <bases>] [--min_match <bases>] [--blat <blat_path>] [--refgene <refGene_path>] [--hg19 <hg19_path>] [--tabix <tabix_path>] [--samtools <samtools_path>] [--prefix <prefix>] [--help]

Description:
This program processes DNA/RNA data to verify and annotate fusion breakpoints in genes. It reads from an input file and a BAM file, aligns sequences using BLAT, and annotates breakpoints with tabix and refGene.

Parameters:
  --bam <bam_file>          : Path to the BAM file for alignment.
  --regions <"regions">     : Regions to check, separated by space (e.g., "chr20:43961871 chr6:117650683").
  --prefix <prefix>         : Prefix for the output files (default: $prefix).
  --extension <bases>       : Number of bases to extend when extracting from the BAM file (default: $extend).
  --min_match <bases>       : Minimum number of matching bases to include in fusion determination (default: $match_min).
  --blat <blat_path>        : Path to the BLAT executable (default: $blat).
  --refgene <refGene_path>  : Path to the refGene database (default: $refgene).
  --hg19 <hg19_path>        : Path to the hg19 reference genome (default: $hg19).
  --tabix <tabix_path>      : Path to the tabix executable (default: $tabix).
  --samtools <samtools_path>: Path to the samtools executable (default: $samtools).
  --out_dir <out_dir>       : Output directory path, default is the current directory.
  --help                    : Display this help message.

Functionality:
1. Data Verification:
   - Reads DNA/RNA data from the input file and verifies the integrity and consistency of the data.
   
2. Fusion Event Annotation:
   - Uses BLAT to align sequences from the BAM file to a reference genome.
   - Annotates fusion breakpoints using tabix and refGene to identify associated genes.

Examples:
  Basic usage:
    perl $0 --bam alignments.bam --regions "chr20:43961871 chr6:117650683" \\
        --extension 100 --min_match 50 \\
        --blat /path/to/blat \\
        --refgene /path/to/refGene.txt \\
        --hg19 /path/to/hg19.fa \\
        --tabix /path/to/tabix \\
        --samtools /path/to/samtools


HELP
    exit;
}

my $start_time = [gettimeofday];

my @breakpoint = split (" ", $regions);
my %default = (
    blat => $blat,
    hg19 => $hg19,
    samtools => $samtools,
    tabix => $tabix,
    refGene => $refgene,
    bam => $bam_file,
    # regions => $regions
);
logger(\%default);
my @cutoff = (
    "extension: $extend",
    "min_match: $match_min",
    "prefix: $prefix",
);
logger(\@cutoff);
for(my $i=0;$i<=$#breakpoint;$i++){
    logger ("index: $i, value: $breakpoint[$i]");
}

die logger("Compressed refGene file or index file not found. Make sure to run bgzip and tabix.\n") if (!-e $refgene||!-e "$refgene.tbi");

my $region;
my @chrom_list;
for (values @breakpoint){
    die logger("Invalid location format. Expected format: chr<chromosome>:<position>, e.g., chr17:57697534\n") unless ($_ =~ /^chr[\w\d]+:\d+$/);
    my ($chrom, $position) = split /:/, $_;
    my $left = $position - $extend;
    my $right = $position + $extend;
    $region .= "$chrom:$left-$right ";
    push @chrom_list, $chrom;
}

logger("regions: $region");
my %seen;
my @chrom_uniq = grep { !$seen{$_}++ } @chrom_list;
# say @chrom_list;
# say @chrom_uniq;

# =pod
logger("Extract potential reads within the query regions from the BAM file.");
open my $target_reads, "-|", "$samtools view $bam_file $region" or die "Could not open pipe: $!";
my $fasta_out = "$out_dir/${prefix}_target.reads.fa";
open TARGET_READS, ">$fasta_out" or die $!; 
my @reads;
while (my $line = <$target_reads>){
    chomp $line;
    my @fields = split "\t", $line;
    my $read_name = $fields[0];
    my $cigar = $fields[5];
    my $sequence = $fields[9];
    # next if(! $cigar =~ /[NS]/);
    next unless($cigar =~ /N|S/);
    # say "$read_name $sequence @reads";
    if(!grep {$_ eq $sequence} @reads){ ## 相同的碱基分布只保留一条
        push @reads, $sequence;
        print TARGET_READS ">$read_name\n$sequence\n";
    }
}
logger("Write to output file $fasta_out");
close TARGET_READS;
close $target_reads;
# =cut


my $chrom_args = join(",", @chrom_uniq);
my $chrom_length = @chrom_uniq;
my %fusions;
logger("Align the target sequences to the reference genome using BLAT");
open my $result_out, ">$out_dir/${prefix}_read.scan.xls" or die $!;
say $result_out "read_name\tread_length\tfusion_name\tbreakpoint1\tbreakpoint2";
my ($breakpoints, $result) = fusion_single($chrom_args);
my %gene_match;
foreach my $pos (@{$breakpoints}){
    my $gene = refgene_search($pos);
    $gene_match{$pos} = $gene;
    # say "$pos $gene";
}
while (my ($name, $info) = each %{$result}){
    my $break1 = $info->{"breakpoint1"};
    my $break2 = $info->{"breakpoint2"};
    my $gene1 = $gene_match{$break1};
    my $gene2 = $gene_match{$break2};
    my $strand1 = $info -> {"strand1"};
    my $strand2 = $info -> {"strand2"};
    my $region1 = $info -> {"region1"};
    my $region2 = $info -> {"region2"};
    my $bc1 = $strand1 eq "+" ? (split "-",$region1)[1] : (split "-",$region1)[0];
    my $bc2 = $strand2 eq "+" ? (split "-",$region2)[0] : (split "-",$region2)[1];
    my $read_length = $info -> {"q_size"};
    my $fusion_name = "$gene1-$gene2";
    say "$name\t$read_length\t$fusion_name\t$break1:$gene1($strand1/$info->{'region1'})\t$break2:$gene2($strand2/$info->{'region2'})";
    # next unless(abs($bc2-$bc1)<11);
    next unless($gene1 ne $gene2);
    
    if (!exists $fusions{$fusion_name}){
        $fusions{$fusion_name} = {
            fusion => $fusion_name
        };
    }
    push @{$fusions{$fusion_name}->{"reads"}}, $name;
    push @{$fusions{$fusion_name}->{"breakpoint1"}}, $break1;
    push @{$fusions{$fusion_name}->{"breakpoint2"}}, $break2;
    
    say $result_out "$name\t$read_length\t$fusion_name\t$break1:$gene1($strand1/$info->{'region1'})\t$break2:$gene2($strand2/$info->{'region2'})"; # if($gene1 ne $gene2);
}
logger("Write to output file $out_dir/${prefix}_read.scan.xls");

logger("Process alignment results and output potential fusion events");
open my $fusion_out, ">$out_dir/${prefix}_fusion.check.xls";
say $fusion_out "fusion_name\tbreakpoint1\tbreakpoint2\treads_count";
foreach my $fusion_ (sort {@{$b->{"reads"}}<=>@{$a->{"reads"}}} values %fusions){
    my $fusion_name = $fusion_ -> {"fusion"};
    my @break1_sort = sort {(split /:/, $a)[1]<=>(split /:/, $b)[1]} @{$fusion_ -> {"breakpoint1"}};
    my @break2_sort = sort {(split /:/, $b)[1]<=>(split /:/, $a)[1]} @{$fusion_ -> {"breakpoint2"}};
    my $read_number = @{$fusion_ -> {"reads"}};
    # logger("$fusion_name @break1_sort");
    say $fusion_out "$fusion_name\t$break1_sort[0]\t$break2_sort[0]\t$read_number";
    logger("$fusion_name\t$break1_sort[0]\t$break2_sort[0]\t$read_number");
}
logger("Write to output file $out_dir/${prefix}_fusion.check.xls");
close $result_out;
close $fusion_out;

my $stop_time = [gettimeofday];
my $elapsed_time = sprintf("%.2f", tv_interval($start_time));
my $elapsed_minutes = sprintf("%.2f", $elapsed_time / 60);
logger("Total elapsed time: $elapsed_minutes minutes(${elapsed_time}s)");


sub fusion_single{
    my $chrom_args = $_[0];
    my @breakpoints;
    my %reads_tmp;
    my %fusions_tmp;
    my %result;
    my $blat_out = "./${prefix}_target.blat.psl";
    my $cmd = "$blat $blat_args $hg19:$chrom_args $fasta_out $blat_out";
    logger($cmd);
    `$cmd`;  ## 上线请取消
    open my $blat_psl, "<$blat_out" or die "open file failed: $!";
    my $count = 0;
    while (<$blat_psl>){
        chomp;
        # next unless($_ =~ /^psL|^\s|^match/);
        next unless(/^\d+/);
        # say;
        my @line = split;
        my $name = $line[9];
        my $t_name = $line[13];
        my $block = $line[17];
        my @block_size = split ",", $line[18];
        my $match = $line[0];
        my $mismatch = $line[1];
        my @t_starts = split ",", $line[-1];
        my $q_size = $line[10];
        my @q_starts = split ",", $line[19];
        my $strand = $line[8];
        $count ++;
        @{$reads_tmp{$count}} = @line;
        # $, = ",";
        # say %reads_tmp;
        # say "@{$reads_tmp{$count}}";
        if ($block == 2){
            next unless($q_size == $match);  ## 适合融合基因在同一条染色体
            $result{$name} = {
                # name => $name,
                strand1 => $strand,
                strand2 => $strand,
                q_size => $q_size,
                breakpoint1 => $t_name.":".$t_starts[0],
                block1 => $block_size[0],
                region1 => region_adjust($strand, $q_starts[0], $block_size[0], $q_size), #($q_starts[0]+1)."-".($q_starts[0]+$block_size[0]),
                breakpoint2 => $t_name.":".$t_starts[1],
                block2 => $block_size[1],
                region2 => region_adjust($strand, $q_starts[1], $block_size[1], $q_size), #($q_starts[1]+1)."-".($q_starts[1]+$block_size[1]),
            };
            push @breakpoints, $t_name.":".$t_starts[0] if(!grep {$_ eq $t_name.":".$t_starts[0]} @breakpoints);
            push @breakpoints, $t_name.":".$t_starts[1] if(!grep {$_ eq $t_name.":".$t_starts[1]} @breakpoints);
        } else{
            next unless($match > $match_min);
            next unless($block == 1);
            # say $name;
            if (!exists $fusions_tmp{$name}){
                $fusions_tmp{$name} = {
                    read_name => $name,
                };
            }
            my $breakpoint = $t_name.":".$t_starts[0];
            push @{$fusions_tmp{$name} -> {"breakList"}}, $breakpoint;
            push @{$fusions_tmp{$name} -> {"rows"}}, $count;
            # say "$name $breakpoint $count";
        }
    }
    foreach my $read_tmp (keys %fusions_tmp){
        my @breakList = @{$fusions_tmp{$read_tmp} -> {"breakList"}};
        my $break_size = @breakList;
        my @rows = @{$fusions_tmp{$read_tmp} -> {"rows"}};
        # say "$read_tmp $break_size";
        next unless($break_size == 2);
        my $breakpoint1 = $breakList[0];
        my $row1 = $rows[0];
        my $breakpoint2 = $breakList[1];
        my $row2 = $rows[1];
        my @line_1 = @{$reads_tmp{$row1}};
        my @line_2 = @{$reads_tmp{$row2}};
        # say $row1;
        # line_1
        my @block_size1 = split ",", $line_1[18];
        my $match1 = $line_1[0];
        my $mismatch1 = $line_1[1];
        my @t_starts1 = split ",", $line_1[-1];
        my @q_starts1 = split ",", $line_1[19];
        my $strand1 = $line_1[8];
        # line_2
        my @block_size2 = split ",", $line_2[18];
        my $match2 = $line_2[0];
        my $mismatch2 = $line_2[1];
        my @t_starts2 = split ",", $line_2[-1];
        my @q_starts2 = split ",", $line_2[19];
        my $strand2 = $line_2[8];
        logger("$read_tmp $breakList[0]:$block_size1[0] $breakList[1]:$block_size2[0]");
        $result{$read_tmp} = {
            strand1 => $strand1,
            strand2 => $strand2,
            q_size => $line_1[10],
            breakpoint1 => $breakpoint1,
            block1 => $block_size1[0],
            region1 => region_adjust($strand1, $q_starts1[0], $block_size1[0], $line_1[10]), #($q_starts1[0]+1)."-".($q_starts1[0]+$block_size1[0]),
            breakpoint2 => $breakpoint2,
            block2 => $block_size2[0],
            region2 => region_adjust($strand2, $q_starts2[0], $block_size2[0], $line_2[10]),
        };
        push @breakpoints, $breakpoint1 if (!grep {$_ eq $breakpoint1} @breakpoints);
        push @breakpoints, $breakpoint2 if (!grep {$_ eq $breakpoint2} @breakpoints);
    }
    close $blat_psl;
    # say "breakpoints: @breakpoints";
    return (\@breakpoints, \%result);
}

sub region_adjust{
    my ($strand, $q_start, $block_size, $q_size) = @_;
    if ($strand eq "+"){
        return ($q_start+1)."-".($q_start+$block_size);
    } elsif($strand eq "-"){
        return ($q_size-$q_start)."-".($q_size-$q_start-$block_size);
    }else{
        die logger("Invalid strand, only '+' and '-' are recognized");
    }
}

sub logger{
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    $year += 1900;
    $mon += 1;
    my $formatted_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec);
    my ($infos) = @_;
    if (ref($infos) eq "HASH"){
        foreach my $key (keys %$infos){
            my $info = $infos -> {$key};
            say "[$formatted_time] - FusionCheck: $key: $info";
        }
    } elsif (ref($infos) eq "ARRAY"){
        foreach my $info (@$infos){
            say "[$formatted_time] - FusionCheck: $info";
        }
    } else{
        say "[$formatted_time] - FusionCheck: $infos";
    }
}

sub position_extend{
    my ($chrom, $position) = split /:/, $_[0];
    my $extend = $_[1];
    my $left = $position - $extend;
    my $right = $position + $extend;
    return "$chrom:$left-$right";
}

sub refgene_search{
    my $pos = $_[0];
    my $pos_extend = position_extend($pos, 1);
    my $cmd = "$tabix $refgene $pos_extend";
    my $gene = ".";
    open my $ref_reader, "-|", $cmd;
    while(<$ref_reader>){
        chomp;
        my @rows = split;
        my $strand = $rows[3];
        $gene = $rows[12];
        last;
    }
    close $ref_reader;
    return $gene;
}

=pod
0      match:    99
1  mis-match:    0
2  rep.match:    0
3        N's:    0
4Q gap count:    1
5Q gap bases:    1
6T gap count:    2
7T gap bases:    4365
8     strand:    +
9     Q name:    V350260186L1C005R02200108900
10     Q size:    100
11    Q start:    0
12      Q end:    100
13     T name:    chr17
14     T size:    81195210  ## 染色体 chr17 全长 
15    T start:    57733325
16      T end:    57737789
17block count:    3
18 blockSizes:    11,50,38,
19    qStarts:    0,12,62,
20    tStarts:    57733325,57733338,57737751,
=cut