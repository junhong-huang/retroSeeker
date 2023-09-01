#!/usr/bin/perl -w

BEGIN {
    unshift( @INC, "/public/home/huangjh/perlBioTools/packages" );
}
use Getopt::Long;
use strict;
use commonTools;
use bioUtils;
use fastTools;
use alignTools;
use samTools;
use annoTools;
use evolutionTools;
use XML::Simple;

$| = 1;


our $queryFile="./mm10.fa";
our $targetFile="./hg38.fa";
my $infoFile  = "";
my $outputDir = "./";
my $help      = 0;
my $debug     = 0;

if ( scalar(@ARGV) < 2 ) {
    usage();
}

my $optLong = GetOptions(
    "query_fasta=s"   => \$queryFile,    # configure file
    "target_fasta=s"   => \$targetFile,    # configure file
    "output=s" => \$outputDir,           # outputDir
    "help|h"   => \$help,
    "debug|d"  => \$debug
) or usage();

&usage() if $help;

sub usage
### usage function
{
    print
"Usage: perl $0 [options] --info <information file> --output <output file>\n";
    print "--info     information file includes genomes for lastzChainNet\n";
    print "--conf     configure file[default=lastzChainNetConfigure.xml]\n";
    print "--output   output dir name[default=axtChainNet]\n";
    print "--help     help information\n";
    print "--debug    run the debug section\n";
    exit;
}

# if ( $infoFile eq "" ) {
#     die " Please set the option --info <information file>\n";
# }

if ( !( -e $outputDir ) ) {
    warn "ceate output dir\n";
    &executeCommand("mkdir $outputDir");
}


my $command_line_example="
perl make_pairwise_alignment_pipeline.pl --query_fasta ./hg38.fa --target_fasta ./mm10.fa --output ./output_data
";

warn "#### two genome->lastz->axt->axtChain->chain->chainNet->net->(syn)";

&doLastzChainNetNew( $outputDir );

sub doLastzChainNetNew {
    my ( $outputDir ) = @_;
    

    my $queryGenomeDir  = "$outputDir/queryGenomeDir";
    (my $queryOrg  = basename($queryFile))=~s/\.fa//;
    my $targetGenomeDir = "$outputDir/targetGenomeDir";
    (my $targetOrg  = basename($targetFile))=~s/\.fa//;

    my $prefix = $targetOrg . "." . $queryOrg;

    my $liftOverDir = $outputDir . "/" . "liftOver";
    if ( !( -e $liftOverDir ) ) {
        warn "ceate output dir\n";
        &executeCommand("mkdir $liftOverDir");
    }

    my $axtOutDir = $outputDir . "/" . "axtChainNet";
    if ( !( -e $axtOutDir ) ) {
        warn "ceate output dir\n";
        &executeCommand("mkdir $axtOutDir");
    }

    my $newOutDir = $axtOutDir . "/" . $targetOrg . "_vs_" . $queryOrg;
    if ( !( -e $newOutDir ) ) {
        warn "ceate output dir\n";
        &executeCommand("mkdir $newOutDir");
    }

    my $netFile =
      &generateAxtChainNet( $targetGenomeDir, $targetOrg,
        $queryGenomeDir, $queryOrg, $newOutDir, $liftOverDir, $prefix );


}

sub generateAxtChainNet {
    my ( $targetGenomeDir, $targetOrg, $queryGenomeDir, $queryOrg, $outputDir, $liftOverDir, $prefix )= @_;

    my $targetChromSizeFile = $targetGenomeDir . "/" . $targetOrg . "/" . $targetOrg . ".chrom.sizes";
    if(1){
        executeCommand("samtools faidx $targetFile");
        executeCommand("cut -f1,2 $targetFile.fai > $targetChromSizeFile");
    }
    my $queryChromSizeFile = $queryGenomeDir . "/" . $queryOrg . "/" . $queryOrg . ".chrom.sizes";
    if(1){
        executeCommand("samtools faidx $queryFile");
        executeCommand("cut -f1,2 $queryFile.fai > $queryChromSizeFile");
    }

    my $outputAxtFile = $outputDir . "/" . $prefix . ".lastz.axt";
    my $commandLine =
        "lastz" . " " 
      . $targetFile . " "
      . $queryFile . " "
      . "--format=axt --ambiguous=iupac ‑‑action:target=multiple --strand=both --allocate:traceback=1.99G" . " > "
      . $outputAxtFile;
    &executeCommand($commandLine);

    my $outputChainFile = $outputDir . "/" . $prefix . ".lastz.axt.chain";
    $commandLine =
        "axtChain" ### download from USCS
      . " -linearGap=loose "
      . $outputAxtFile
      . " -faT "
      . $targetFile
      . " -faQ "
      . $queryFile . " "
      . $outputChainFile;
    &executeCommand($commandLine);

    my $outputTargetNetFile = $outputDir . "/" . $prefix . ".lastz.axt.chain.target.net";
    my $outputQueryNetFile = $outputDir . "/" . $prefix . ".lastz.axt.chain.query.net";
    $commandLine =
        'chainNet' . " " ### download from USCS
      . $outputChainFile . " "
      . $targetChromSizeFile . " "
      . $queryChromSizeFile . " "
      . $outputTargetNetFile . " "
      . $outputQueryNetFile;
    &executeCommand($commandLine);

    my $outputTargetNetSyntenicFile = $outputDir . "/" . $prefix . ".lastz.axt.chain.target.syntenic.net";
    $commandLine =
        'netSyntenic' . " " ### download from USCS
      . $outputTargetNetFile . " "
      . $outputTargetNetSyntenicFile;
    &executeCommand($commandLine);

    my $outputQueryNetSyntenicFile = $outputDir . "/" . $prefix . ".lastz.axt.chain.query.syntenic.net";
    $commandLine =
        'netSyntenic' . " " ### download from USCS
      . $outputQueryNetFile . " "
      . $outputQueryNetSyntenicFile;
    &executeCommand($commandLine);

    my $outputSubsetChainFile = $outputDir . "/" . $prefix . ".lastz.axt.subset.net.chain";
    $commandLine =
        'netChainSubset' . " " ### download from USCS
      . $outputTargetNetFile . " "
      . $outputChainFile . " "
      . $outputSubsetChainFile;
    &executeCommand($commandLine);

    my $outputOverChainFile = $liftOverDir . "/" . $prefix . ".over.chain";
    $commandLine =
        'chainStitchId' . " " ### download from USCS
      . $outputSubsetChainFile . " "
      . $outputOverChainFile;
    &executeCommand($commandLine);

    return $outputOverChainFile;
}
