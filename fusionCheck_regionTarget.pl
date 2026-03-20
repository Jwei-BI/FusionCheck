#!/usr/bin/perl
use strict;
use Cwd;
use File::Spec;
use Getopt::Long;
use feature qw(say);
use List::Util qw(sum);
# use Digest::MD5 qw(md5_hex);
use Time::HiRes qw(gettimeofday tv_interval);
# use JSON;

# my $json = JSON->new->utf8->pretty;

# Default parameter values
my $blat        = "/software/blat";
my $blat_args   = "-fastMap -minMatch=2 -maxGap=1 -tileSize=11 -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -t=dna -q=dna -out=psl";
# my $blat_args   = "-fastMap -out=psl";
my $hg19        = "/ref/hg19.2bit";
my $samtools    = "/software/samtools/samtools-1.9/samtools";
my $tabix       = "/software/samtools/samtools-1.9/htslib-1.9/tabix";
my $refgene     = "/ref/refGene.txt.sort.gz";
my $rmsk        = "/ref/hg19.rmsk.sort.gz";
my $extend      = 100;
my $match_min   = 20;
my $out_dir     = getcwd();
my $extend_pair = 1000;
my $break_region = 10000;
my $regions;
my $bam_file;
my $bed_file;
my $prefix = "fusionCheck";
my $fast;
my $filter;
my $help;

# Parse command line parameters
GetOptions(
    'bam=s'         => \$bam_file,
    'regions=s'     => \$regions,
    'bed=s'         => \$bed_file,
    'extension=i'   => \$extend,
    'min_match=i'   => \$match_min,
    'blat=s'        => \$blat,
    'samtools=s'    => \$samtools,
    'tabix=s'       => \$tabix,
    'hg19=s'        => \$hg19,
    'refgene=s'     => \$refgene,
    'rmsk=s'        => \$rmsk,
    'prefix=s'      => \$prefix,
    'out_dir=s'     => \$out_dir,
    'fast'          => \$fast,
    'filter'        => \$filter,
    'help'          => \$help,
) or die "Incorrect usage!\n";

# Display help information
if ($help || !$bam_file || (!$regions and !$bed_file)) {
    print <<HELP;
Usage: $0 --bam <bam_file> [--regions <"region1 region2"> or --bed <bed_file>] [--filter --fast]
        [--extension <bases>] [--min_match <bases>] [--blat <blat_path>] [--refgene <refGene_path>] [--rmsk <rmsk_path>]
        [--hg19 <hg19_path>] [--tabix <tabix_path>] [--samtools <samtools_path>] [--prefix <prefix>] [--help]

Description:
This program processes DNA/RNA data to verify and annotate fusion breakpoints in genes.
It reads from an input file and a BAM file, aligns sequences using BLAT, and annotates breakpoints with tabix and refGene.

Parameters:
  --bam <bam_file>          : Path to the BAM file for alignment.
  --regions <"regions">     : Regions to check, separated by space (e.g., "chr20:43961871 chr6:117650683").
  --bed <bed_file>          : BED file specifying the regions to scan for fusion events.
  --prefix <prefix>         : Prefix for the output files (default: $prefix).
  --extension <bases>       : Number of bases to extend when extracting from the BAM file (default: $extend).
  --min_match <bases>       : Minimum number of matching bases to include in fusion determination (default: $match_min).
  --blat <blat_path>        : Path to the BLAT executable (default: $blat).
  --refgene <refGene_path>  : Path to the refGene database (default: $refgene).
  --rmsk <rmsk_path>        : Path to the rmsk database (default: $rmsk).
  --hg19 <hg19_path>        : Path to the hg19 reference genome (default: $hg19).
  --tabix <tabix_path>      : Path to the tabix executable (default: $tabix).
  --samtools <samtools_path>: Path to the samtools executable (default: $samtools).
  --out_dir <out_dir>       : Output directory path, default is the current directory.
  --fast                    : Fast model(.2bit).
  --filter                  : Enable fusion filter module and write to .fusion.filter.xls.";
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
  If exists bed:
    perl $0 --bam alignments.bam --bed <bed_file> \\
        --extension 100 --min_match 50 \\
        --blat /path/to/blat \\
        --refgene /path/to/refGene.txt \\
        --hg19 /path/to/hg19.fa \\
        --tabix /path/to/tabix \\
        --samtools /path/to/samtools

author = Ma Junwei
  date = 2024-05-30
  mail = 2254170457\@qq.com

HELP
    exit;
}

die logger("Compressed refGene file or index file not found. Make sure to run bgzip and tabix.\n") if (!-e $refgene||!-e "$refgene.tbi");
die logger("Compressed rmsk file or index file not found. Make sure to run bgzip and tabix.\n") if (!-e $rmsk||!-e "$rmsk.tbi");
my %default = (
    blat => $blat,
    hg19 => $hg19,
    samtools => $samtools,
    tabix => $tabix,
    refGene => $refgene,
    bam => $bam_file,
);
logger(\%default);
my @cutoff = (
    "extension: $extend",
    "min_match: $match_min",
    "prefix: $prefix",
);
logger(\@cutoff);

my $start_time = [gettimeofday];
my $region;
my $region_blat;
my $target_reads;
# my $fasta_out = "$out_dir/${prefix}_target.reads.fa";
my @regions_list;
# my @regions_blat_list;
logger("Extract potential reads within the query regions from the BAM file.");
if (defined($regions)){
    my @breakpoint = split (" ", $regions);

    for (values @breakpoint){
        die logger("Invalid location format. Expected format: chr<chromosome>:<position>, e.g., chr17:57697534\n") unless ($_ =~ /^chr[\w\d]+:\d+$/);
        my ($chrom, $position) = split /:/, $_;
        my $left = $position - $extend;
        my $right = $position + $extend;
        $region .= "$chrom:$left-$right ";
        $region_blat .= "$chrom:$left-$right,";
    }

    # logger("regions: $region");
    # logger("regions_blat: $region_blat");
    push @regions_list, $region;
    # push @regions_blat_list, $region_blat;
    # open $target_reads, "-|", "$samtools view $bam_file $region" or die "Could not open pipe: $!";
} elsif(defined($bed_file)){
    logger("bed_file: $bed_file");
    # open $target_reads, "-|", "$samtools view $bam_file -L $bed_file -q 1" or die "Could not open pipe: $!";
    open my $bed_reader, "<$bed_file" or die "open file failed: $!";
    while(my $line = <$bed_reader>) {
        chomp $line;
        my @fields = split /\t/, $line;
        my ($chrom1, $start1, $stop1, $gene1, $strand1,
            $chrom2, $start2, $stop2, $gene2, $strand2) = @fields;
        my $left = $chrom1.":".($start1-$extend)."-".($stop1+$extend);
        my $right = $chrom2.":".($start2-$extend)."-".($stop2+$extend);
        push @regions_list, "$left $right";
        # push @regions_blat_list, "$left,$right";
    }
}

# open TARGET_READS, ">$fasta_out" or die $!; 
# my @reads;
# while (my $line = <$target_reads>){
#     chomp $line;
#     my @fields = split "\t", $line;
#     my $read_name = $fields[0];
#     my $cigar = $fields[5];
#     my $sequence = $fields[9];
#     my $flag = $fields[1];
#     my $pair;
#     if ($flag & 0x40){
#         $pair = 1;
#     } elsif($flag & 0x80){
#         $pair = 2;
#     }else{
#         $pair = ".";
#     }
#     next unless($cigar =~ /N|S/);

#     if(!grep {$_ eq $sequence} @reads){
#         push @reads, $sequence;
#         print TARGET_READS ">$read_name/$pair\n$sequence\n";
#     }
# }
# logger("Writing to output file $fasta_out");
# close TARGET_READS;
# close $target_reads;

my %filter_hash;

my %gene_match;
my %gene_strand;
my %reads_align;
my $region_count;
foreach my $region (@regions_list){
    $region_count ++;
    # $tmp_dir =~ s/ /_/g;
    # my $hash_ = md5_hex($tmp_dir);
    my $tmp_hex = "tmp_$region_count";
    process_cmd("mkdir -p $out_dir/tmp/$tmp_hex");
    my $tmp = "$out_dir/tmp/$tmp_hex";
    my %fusions;
    logger("Align the target sequences to the reference genome using BLAT");
    open my $result_out, ">$tmp/${prefix}_read.scan.xls" or die $!;
    say $result_out "read_name\tread_length\tfusion_name\tbreakpoint1\tbreakpoint2\tscore";
    # my ($breakpoints, $result) = fusion_single($chrom_args);
    my $region_blat = $region;
    $region_blat =~ s/ /,/g;
    logger("#$region_count regions: $region");
    (!$fast && $region_count==1) ? logger("#$region_count regions_blat: all")
                                 : logger("#$region_count regions_blat: $region_blat");
    my ($breakpoints, $result) = fusion_single($region, $region_blat);

    foreach my $pos (@{$breakpoints}){
        if (!exists $gene_match{$pos}){
            my ($gene, $strand) = refgene_search($pos);
            $gene_match{$pos} = $gene;
            $gene_strand{$gene} = $strand;
        }
    }
    while (my ($name, $info) = each %{$result}){
        my $break1 = $info->{"breakpoint1"};
        my ($chrom1, $bp1) = split /:/, $break1;
        my $break2 = $info->{"breakpoint2"};
        my ($chrom2, $bp2) = split /:/, $break2;
        my $gene1  = $gene_match{$break1};
        my $gene2  = $gene_match{$break2};
        my $strand1 = $info -> {"strand1"};
        my $strand2 = $info -> {"strand2"};
        my $region1 = $info -> {"region1"};
        my $region2 = $info -> {"region2"};
        my $bc1 = $strand1 eq "+" ? (split "-",$region1)[1] : (split "-",$region1)[0];
        my $bc2 = $strand2 eq "+" ? (split "-",$region2)[0] : (split "-",$region2)[1];
        my $read_length = $info -> {"q_size"};
        my $score = $info -> {"score"};
        next unless(abs($bc2-$bc1)<11);
        next if(abs($bp1-$bp2)<$break_region&&$chrom1 eq $chrom2);
        next unless($gene1 ne $gene2);
        my $direction_left = $strand1 eq "+" ? "reversed" : "forward";
        my $direction_right = $strand2 eq "+" ? "forward" : "reversed";
        my $fusion_name;
        my $fusion_gene;
        if ($gene1 ge $gene2){
            $fusion_gene = "$gene1-$direction_left"."::"."$gene2-$direction_right";
            $fusion_name = "$break1-$direction_left"."::"."$break2-$direction_right";
        } else{
            $fusion_gene = "$gene2-$direction_right"."::"."$gene1-$direction_left";
            $fusion_name = "$break2-$direction_right"."::"."$break1-$direction_left";
            $break1 = $info->{"breakpoint2"};
            $break2 = $info->{"breakpoint1"};
            $gene1  = $gene_match{$break1};
            $gene2  = $gene_match{$break2};
            $strand1 = $info -> {"strand2"};
            $strand2 = $info -> {"strand1"};
            $region1 = $info -> {"region2"};
            $region2 = $info -> {"region1"};
        }
        # my $fusion_name = "$gene1-$direction_left"."::"."$gene2-$direction_right";
        if (!exists $fusions{$fusion_name}){
            $fusions{$fusion_name} = {
                fusion => $fusion_gene,
                break_name => $fusion_name
            };
        }
        push @{$fusions{$fusion_name}->{"reads"}}, $name;
        push @{$fusions{$fusion_name}->{"breakpoint1"}}, $break1;
        push @{$fusions{$fusion_name}->{"breakpoint2"}}, $break2;
        push @{$fusions{$fusion_name}->{"scores"}}, $score;

        say $result_out "$name\t$read_length\t$fusion_gene\t"
                    . "$break1:$gene1($strand1/$region1)\t"
                    . "$break2:$gene2($strand2/$region2)\t"
                    . "$score";
        $reads_align{$name} = "$name\t$read_length\t$fusion_gene\t"
                            . "$break1:$gene1($strand1/$region1)\t"
                            . "$break2:$gene2($strand2/$region2)\t"
                            . "$score" if (!exists $reads_align{$name});
    }
    # print $json->encode(\%fusions);
    logger("Writing to output file $tmp/${prefix}_read.scan.xls");

    logger("Process alignment results and output potential fusion events");
    open my $fusion_out, ">$tmp/${prefix}.fusion.xls";

    say $fusion_out "fusion_name\tbreakpoint1\tbreakpoint2\treads_count\tspanning_count\tscore_mean";
    open my $spanning_out, ">$tmp/${prefix}.spanning.xls";
    say $spanning_out "fusion_name\tspanning_reads";
    foreach my $fusion_ (sort {@{$b->{"reads"}}<=>@{$a->{"reads"}}} values %fusions){
        my $fusion_name = $fusion_ -> {"fusion"};
        my $break_name = $fusion_ -> {"break_name"};
        my @break1_sort = sort {(split /:/, $b)[1]<=>(split /:/, $a)[1]} @{$fusion_ -> {"breakpoint1"}};
        my @break2_sort = sort {(split /:/, $a)[1]<=>(split /:/, $b)[1]} @{$fusion_ -> {"breakpoint2"}};
        my @reads = @{$fusion_ -> {"reads"}};
        my $read_number = @reads;
        my @scores      = @{$fusion_ -> {"scores"}};
        my $score_sum   = sum(@scores);
        my $score_mean  = $read_number ? sprintf "%.2f", $score_sum/$read_number : 0;
        my %seen1;
        my %seen2;
        my @chroms_1 = grep {!$seen1{$_}++} (map{(split /:/, $_)[0]} @break1_sort);
        my @chroms_2 = grep {!$seen2{$_}++} (map{(split /:/, $_)[0]} @break2_sort);
        my $chrom_num1 = @chroms_1;
        my $chrom_num2 = @chroms_2;
        next unless($chrom_num1==1 && $chrom_num2==1);
        my ($block1, $block2) = split /::/, $fusion_name;
        my ($gene1, $direction1) = $block1 =~ /(.*)-(.*)/;
        my ($gene2, $direction2) = $block2 =~ /(.*)-(.*)/;
        my @spanning = discordant_pairs($break1_sort[0], $break2_sort[0], $direction1, $direction2);
        my $spanning_num = @spanning;
        my $spanning_reads = $spanning_num>0 ? join(",", @spanning) : ".";
        say $spanning_out "$fusion_name\t$spanning_reads";
        say $fusion_out "$fusion_name\t$break1_sort[0]\t$break2_sort[0]\t$read_number\t$spanning_num\t$score_mean";
        if ($filter){
            next if($break1_sort[0] =~ /_/);
            next if($break2_sort[0] =~ /_/);
            # next unless($read_number>2 && $spanning_num>0);
            # next unless($read_number>2);
            # say $filter_out "$fusion_name\t$break1_sort[0]\t$break2_sort[0]\t$read_number\t$spanning_num\t$score_mean";
            logger("$fusion_name\t$break1_sort[0]\t$break2_sort[0]\t$read_number\t$spanning_num\t$score_mean");
            if (!exists $filter_hash{$break_name}){
                $filter_hash{$break_name} = {
                    fusion => $fusion_name,
                    break_name => $break_name,
                    reads => \@reads,
                    break1 => $break1_sort[0],
                    break2 => $break2_sort[0],
                    spanning_num => $spanning_num,
                    scores => \@scores,
                    spanning_reads => \@spanning
                };
            } else {
                # if (!grep {$_ eq $breakpoint1} @breakpoints)
                foreach my $i (0..$#reads){
                    if (!grep {$_ eq $reads[$i]} @{$filter_hash{$break_name} -> {"reads"}}){
                        push (@{$filter_hash{$break_name} -> {"reads"}}, $reads[$i]);
                        push (@{$filter_hash{$break_name} -> {"scores"}}, $scores[$i]);
                    }
                }
            }
        }
    }
    logger("Writing to output file $tmp/${prefix}.fusion.xls");
    close $result_out;
    close $fusion_out;
    close $spanning_out;
    # close $filter_out if ($filter);

}

my $filter_out;
if ($filter){
    open $filter_out, ">$out_dir/${prefix}.fusion.filter.xls";
    say $filter_out "fusion_name\tbreakpoint1\trmsk1\tbreakpoint2\trmsk2\treads_count\tspanning_count\tscore_mean";
}
open RESULT_OUT, ">$out_dir/${prefix}.fusion.all.xls";
say RESULT_OUT "fusion_name\tbreakpoint1\trmsk1\tbreakpoint2\trmsk2\treads_count\tspanning_count\tscore_mean";
open READ_SCAN, ">$out_dir/${prefix}.read.scan.xls";
say READ_SCAN "read_name\tread_length\tfusion_name\tbreakpoint1\tbreakpoint2\tscore";
open SPANNING_OUT, ">$out_dir/${prefix}.spanning.all.xls";
say SPANNING_OUT "fusion_name\tspanning_count\tspanning_reads";
foreach my $fs (sort {@{$b->{"reads"}}<=>@{$a->{"reads"}}} values %filter_hash){
    my $fusion = $fs -> {"fusion"};
    my $break_name = $fs -> {"break_name"};
    my $break1 = $fs -> {"break1"};
    my $break2 = $fs -> {"break2"};
    my @reads  = @{$fs -> {"reads"}};
    my $read_number = @reads;
    my $spanning_num = $fs -> {"spanning_num"};
    my @scores = @{$fs -> {"scores"}};
    my $score_sum   = sum(@scores);
    my $score_mean  = $read_number ? sprintf "%.2f", $score_sum/$read_number : 0;
    my %seen;
    my @reads_uniq = grep {!$seen{$_}++} @reads;
    my $read_uniq = @reads_uniq;
    my @spanning = @{$fs -> {"spanning_reads"}};
    my $spanning_count = @spanning;
    my $spanning_reads = $spanning_num>0 ? join(",", @spanning) : ".";
    say SPANNING_OUT "$fusion\t$spanning_count\t$spanning_reads";
    # say "$fusion\t$break1\t$break2\t$read_uniq\t$spanning_num\t$score_mean";
    foreach my $read_name (@reads_uniq){
        say READ_SCAN $reads_align{$read_name};
    }
    my ($rmsk1, $rmsk2) = rmsk_annot($break_name);
    say RESULT_OUT "$fusion\t$break1\t$rmsk1\t$break2\t$rmsk2\t$read_uniq\t$spanning_num\t$score_mean";
    if ($filter){
        next unless($read_uniq>2);
        next unless(($rmsk1 eq "." && $rmsk2 eq ".")
                    || ($spanning_num>0)
                    || ($read_uniq>5 && $score_mean>0.98));
        say $filter_out "$fusion\t$break1\t$rmsk1\t$break2\t$rmsk2\t$read_uniq\t$spanning_num\t$score_mean";
    }
}

close READ_SCAN;
close $filter_out if ($filter);

my $stop_time = [gettimeofday];
my $elapsed_time = sprintf("%.2f", tv_interval($start_time));
my $elapsed_minutes = sprintf("%.2f", $elapsed_time / 60);
logger("Total elapsed time: $elapsed_minutes minutes(${elapsed_time}s)");
my $mem_usage = get_memory_usage();
my $mem_g = sprintf "%.2f", $mem_usage/1024;
logger("Memory usage: $mem_usage KB/$mem_g MB");

sub rmsk_annot {
    my ($break_name) = @_;
    # = $block1 =~ /(.*)-(.*)/; /(chr*):([0-9]+)-({forward,reversed})::(chr*):([0-9]+)-({forward,reversed})/
    my ($chrom1,$break1,$direction1,$chrom2,$break2,$direction2) = $break_name =~ /(chr\d+):(\d+)-(forward|reversed)::(chr\d+):(\d+)-(forward|reversed)/;
    # say "$chrom1,$break1,$direction1,$chrom2,$break2,$direction2";
    # my $rmsk_left = ".";
    # my $rmsk_right = ".";
    my $rmsk_left = breakpoint_extend($chrom1, $break1, $direction1);
    my $rmsk_right = breakpoint_extend($chrom2, $break2, $direction2);
    return ($rmsk_left, $rmsk_right)
}

sub breakpoint_extend {
    my ($chrom, $breakpoint, $direction) = @_;
    my $key_ = "$chrom:$breakpoint";
    my $gene_ = $gene_match{$key_};
    my $strand = $gene_strand{$gene_};
    my $region_extend = "$chrom:$breakpoint-".($breakpoint + $extend);
    # my $ex_ = $breakpoint + $extend;
    if ($direction eq "forward"){
        # $ex_ = $breakpoint - $extend;
        $region_extend = "$chrom:".($breakpoint-$extend)."-$breakpoint";
    }
    # return $region_extend;
    my @annot;
    open RMSK_READER, "-|", "$tabix $rmsk $region_extend";
    while (<RMSK_READER>){
        chomp;
        my @line = split;
        my ($chrom_q, $start_q, $stop_q) = ($line[4], $line[5], $line[6]);
        my ($match, $class) = ($line[9], $line[10]);
        push @annot, "$chrom_q:$start_q-$stop_q,$match,$class";
    }
    return (scalar @annot)>0 ? join(";", @annot) : ".";
}

sub discordant_pairs {
    my ($breakpoint1, $breakpoint2, $direction1, $direction2) = @_;
    my @reads_1 = reads_get($breakpoint1, $direction1);
    my @reads_2 = reads_get($breakpoint2, $direction2);
    my %check;
    $check{$_} ++ for @reads_1;
    my @spanning = grep { $check{$_} } @reads_2;
    return @spanning;
}

sub reads_get {
    my ($break, $direction) = @_;
    my ($chrom, $point) = split /:/, $break;
    my $left = $point - $extend_pair;
    my $right =  $point + $extend_pair;
    my @reads;
    open my $bam_read, "-|", "$samtools view $bam_file $chrom:$left-$right" or die "Could not open pipe: $!";
    while (my $line = <$bam_read>){
        chomp $line;
        my @fields = split "\t", $line;
        my $read_name = $fields[0];
        my $tar_chrom = $fields[2];
        my $tar_point = $fields[3];
        my $qual      = $fields[4];
        my $cigar     = $fields[5];
        my $flag      = $fields[1];
        my $seq       = $fields[9];
        my $seq_length = length($seq);
        next if($flag & 0x8);
        next if($cigar =~ /N|S|H/);
        next if($qual<20);
        if($tar_chrom eq $chrom){
            # if ($tar_point > $point || $tar_point+$seq_length < $point){
            #     push @reads, $read_name;
            # }
            if ($direction eq "forward" && $tar_point > $point){
                push @reads, $read_name;
            } elsif($direction eq "reversed" && $tar_point+$seq_length < $point){
                push @reads, $read_name;
            }
        }
    }
    return @reads;

}

sub fusion_single {
    my ($region, $region_blat) = @_;
    my @breakpoints;
    my %reads_tmp;
    my %fusions_tmp;
    my %result;
    # my $blat_out = "./${prefix}_target.blat.psl";
    my $cmd;

    if (length($region)>0 && ($hg19 =~ /\.2bit$/) && $fast){
        open BLAT_PSL, "-|", "$samtools view $bam_file $region"
                       ."| awk -F'\t' -vOFS='\n' '\$6~/S|N/{pair=0;if(and(\$2,64)){pair=1}else if(and(\$2,128)){pair=2};print \">\"\$1\"/\"pair,\$10}'"
                       ."| $blat $blat_args $hg19:$region_blat stdin stdout";
        # say "-|", "$samtools view $bam_file $region"
        #                ."| awk -F'\\t' -vOFS='\\n' '\$6~/S|N/{pair=0;if(and(\$2,64)){pair=1}else if(and(\$2,128)){pair=2};print \">\"\$1\"/\"pair,\$10}'"
        #                ."| $blat $blat_args $hg19:$region_blat stdin stdout";
    }else{
        open BLAT_PSL, "-|", "$samtools view $bam_file $region"
                       ."| awk -F'\t' -vOFS='\n' '\$6~/S|N/{pair=0;if(and(\$2,64)){pair=1}else if(and(\$2,128)){pair=2};print \">\"\$1\"/\"pair,\$10}'"
                       ."| $blat $blat_args $hg19 stdin stdout";
    }
    # process_cmd($cmd);
    # die logger("failed: $cmd") if ($?!=0);
    # open my $blat_psl, "<$blat_out" or die "open file failed: $!";
    # open my $blat_psl, "-|", $cmd or die logger("failed: $cmd");
    my $count = 0;
    while (<BLAT_PSL>){
        chomp;
        next unless(/^\d+/);
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
        my ($chrom, $start, $stop) = split /:|-/, $t_name;
        @{$reads_tmp{$count}} = @line;

        if ($block == 2){
            next unless($q_size == $match);  ## 适合融合基因在同一条染色体
            my $breakpoint1 = $strand eq "-"
                                ? $chrom.":".($start+$t_starts[0]+1)
                                : $chrom.":".($start+$t_starts[0]+$block_size[0]);
            my $breakpoint2 = $strand eq "+"
                                ? $chrom.":".($start+$t_starts[1]+1)
                                : $chrom.":".($start+$t_starts[1]+$block_size[1]);
            my $gap_base    = $line[5];
            my $N_base      = $line[3];
            my $score_tmp   = ($block_size[0]+$block_size[1]) - ($mismatch+$gap_base+$N_base);
            my $ratio       = sprintf "%.2f", $score_tmp/$q_size;
            $result{$name} = {
                score => $ratio,
                strand1 => $strand,
                strand2 => $strand,
                q_size => $q_size,
                breakpoint1 => $breakpoint1,
                block1 => $block_size[0],
                region1 => region_adjust($strand, $q_starts[0], $block_size[0], $q_size),
                breakpoint2 => $breakpoint2,
                block2 => $block_size[1],
                region2 => region_adjust($strand, $q_starts[1], $block_size[1], $q_size),
            };
            push @breakpoints, $breakpoint1 if(!grep {$_ eq $breakpoint1} @breakpoints);
            push @breakpoints, $breakpoint2 if(!grep {$_ eq $breakpoint2} @breakpoints);
        } else{
            next unless($match > $match_min);
            next unless($block == 1);
            if (!exists $fusions_tmp{$name}){
                $fusions_tmp{$name} = {
                    read_name => $name,
                };
            }
            my $breakpoint = $chrom.":".($start+$t_starts[0]);
            push @{$fusions_tmp{$name} -> {"breakList"}}, $breakpoint;
            push @{$fusions_tmp{$name} -> {"rows"}}, $count;
        }
    }
    foreach my $read_tmp (keys %fusions_tmp){
        my @breakList = @{$fusions_tmp{$read_tmp} -> {"breakList"}};
        my $break_size = @breakList;
        my @rows = @{$fusions_tmp{$read_tmp} -> {"rows"}};
        ## 针对比对情况复杂的事件进行处理
        my ($fs_tmp, $score_tmp) = align_complex(\@rows, \%reads_tmp);
        next unless($fs_tmp);
        my ($row1, $row2) = split /::/, $fs_tmp;
        my @line_1 = @{$reads_tmp{$row1}};
        my @line_2 = @{$reads_tmp{$row2}};
        # line_1
        my @block_size1 = split ",", $line_1[18];
        my $match1 = $line_1[0];
        my $mismatch1 = $line_1[1];
        my @t_starts1 = split ",", $line_1[-1];
        my @q_starts1 = split ",", $line_1[19];
        my $strand1 = $line_1[8];
        my $t_name1 = $line_1[13];
        my ($chrom1, $start1, $end1) = split /:|-/, $t_name1;
        my $stop1 = $line_1[16];
        my $breakpoint1 = $strand1 eq "-" ? "$chrom1:".($start1+$t_starts1[0]+1) : $chrom1.":".($start1+$stop1);
        # line_2
        my @block_size2 = split ",", $line_2[18];
        my $match2 = $line_2[0];
        my $mismatch2 = $line_2[1];
        my @t_starts2 = split ",", $line_2[-1];
        my @q_starts2 = split ",", $line_2[19];
        my $t_name2 = $line_2[13];
        my ($chrom2, $start2, $end2) = split /:|-/, $t_name2;
        my $strand2 = $line_2[8];
        my $stop2 = $line_2[16];
        my $breakpoint2 = $strand2 eq "+" ? "$chrom2:".($start2+$t_starts2[0]+1) : $chrom2.":".($start2+$stop2);
        $result{$read_tmp} = {
            score => $score_tmp,
            strand1 => $strand1,
            strand2 => $strand2,
            q_size => $line_1[10],
            breakpoint1 => $breakpoint1,
            block1 => $block_size1[0],
            region1 => region_adjust($strand1, $q_starts1[0], $block_size1[0], $line_1[10]),
            breakpoint2 => $breakpoint2,
            block2 => $block_size2[0],
            region2 => region_adjust($strand2, $q_starts2[0], $block_size2[0], $line_2[10]),
        };
        push @breakpoints, $breakpoint1 if (!grep {$_ eq $breakpoint1} @breakpoints);
        push @breakpoints, $breakpoint2 if (!grep {$_ eq $breakpoint2} @breakpoints);
    }
    close BLAT_PSL;
    # say "@breakpoints";

    # print $json->encode(\%result);
    return (\@breakpoints, \%result);
}

sub process_cmd {
    my ($cmd) = @_;
    logger("CMD: $cmd");
    open my $return, "-|", $cmd or die "Error, CMD: $cmd failed; $!";
    while (my $line = <$return>){
        chomp $line;
        logger($line);
    }
    return
}

sub align_complex {
    my ($rows, $reads) = @_;
    my @sort_rows = sort {
        @{$reads -> {$b}}[0] <=> @{$reads -> {$a}}[0]
    } @{$rows};
    my %seen;
    my @cross = map{
        my $x = $_;
        map {
            my $y = $_;
            unless ($x eq $y || $seen{"$y,$x"}){
                $seen{"$x,$y"} = 1;
                my $x_start = @{$reads->{$x}}[11];
                my $y_start = @{$reads->{$y}}[11];
                ($x_start < $y_start) ? $x."::".$y : $y."::".$x;
            }else {
                ()
            }
        } (scalar @sort_rows)>5 ? @sort_rows[0..5] : @sort_rows
    } (scalar @sort_rows)>5 ? @sort_rows[0..5] : @sort_rows;
    my $len = @cross;
    my %candidate;
    foreach my $tmp (@cross){
        my ($row1, $row2) = split /::/, $tmp;
        my @read1 = @{$reads -> {$row1}};
        my @read2 = @{$reads -> {$row2}};
        my $tmp_check = region_judge(\@read1, \@read2);
        if (defined($tmp_check)){
            $candidate{$tmp} = $tmp_check;
        }
    }
    return () unless(%candidate);
    my @keys_candidate = keys %candidate;
    my @keys_sort = sort {
        $candidate{$b} <=> $candidate{$a}
    } @keys_candidate;

    my $key_max = @keys_sort[0];
    my $score_max = $candidate{$key_max};
    return ($key_max, $score_max);

}

sub region_judge {
    my ($read1_ref, $read2_ref) = @_;
    my @read1 = @$read1_ref;
    my @read2 = @$read2_ref;
    my $match1 = $read1[0];
    my $mismatch1 = $read1[1];
    my $strand1 = $read1[8];
    my $block1 = $read1[18];
    $block1 =~ s/,$//;
    my $qStart1 = $read1[11];
    my $qEnd1   = $read1[12];
    my $gap_base1 = $read1[5];
    my $N_base1 = $read1[3];
    # //
    my $match2 = $read2[0];
    my $mismatch2 = $read2[1];
    my $strand2 = $read2[8];
    my $block2 = $read2[18];
    $block2 =~ s/,$//;
    my $qStart2 = $read2[11];
    my $qEnd2   = $read2[12];
    my $gap_base2 = $read2[5];
    my $N_base2 = $read2[3];
    my $qSize   = $read1[10];
    my $gap_base = $gap_base1+$gap_base2;
    my $N_base = $N_base1+$N_base2;
    my $match_check = $qEnd1 > $qStart2
                        ? ($block1+$block2) - ($mismatch1+$mismatch2+$gap_base+$N_base) - 2*($qEnd1-$qStart2)
                        : ($block1+$block2) - ($mismatch1+$mismatch2+$gap_base+$N_base);

    my $ratio  = sprintf "%.2f", $match_check/$qSize;
    return $ratio > 0.8 ? $ratio : ()
}

sub region_adjust {
    my ($strand, $q_start, $block_size, $q_size) = @_;
    if ($strand eq "+"){
        return ($q_start+1)."-".($q_start+$block_size);
    } elsif($strand eq "-"){
        return ($q_size-$q_start)."-".($q_size-$q_start-$block_size+1);
    }else{
        die logger("Invalid strand, only '+' and '-' are recognized");
    }
}

sub get_memory_usage {
    open(my $fh, '<', '/proc/self/status') or die "Cannot open /proc/self/status: $!";
    my $memory_usage;
    while (<$fh>) {
        if (/^VmRSS:\s+(\d+)\s+kB/) {
            $memory_usage = $1;
            last;
        }
    }
    close($fh);
    return $memory_usage;
}

sub logger {
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    $year += 1900;
    $mon += 1;
    my $formatted_time = sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                            $year, $mon, $mday, $hour, $min, $sec);
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

sub position_extend {
    my ($chrom, $position) = split /:/, $_[0];
    my $extend = $_[1];
    my $left = $position - $extend;
    my $right = $position + $extend;
    return "$chrom:$left-$right";
}

sub refgene_search {
    my $pos = $_[0];
    my $pos_extend = position_extend($pos, 1);
    my $cmd = "$tabix $refgene $pos_extend";
    my $gene = ".";
    my $strand = ".";
    my $gene_strand = "./.";
    open my $ref_reader, "-|", $cmd;
    while(<$ref_reader>){
        chomp;
        my @rows = split;
        $strand = $rows[3];
        $gene = $rows[12];
        my $cds_start_stat = $rows[13];
        my $cds_end_stat = $rows[14];
        next unless($cds_start_stat eq "cmpl" && $cds_end_stat eq "cmpl");
        $gene_strand = "$gene/$strand";
        last;
    }
    close $ref_reader;
    return ($gene, $strand);
}

