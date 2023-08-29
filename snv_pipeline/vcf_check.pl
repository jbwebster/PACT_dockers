#!/usr/bin/perl

use strict;
use warnings;

use File::Copy qw(copy);
use IPC::Run qw(run);


if(@ARGV < 3) {
    print STDERR "Usage: <VCF> <output destination if VCF empty> <command if VCF non-empty> [command args...]\n";
    exit -1;
}

main(@ARGV);
exit 0;

sub main {
    my ($vcf_to_check, $where_to_copy_if_empty, @commandline_if_nonempty) = @_;

    if (vcf_is_empty($vcf_to_check)) {
        copy_file($vcf_to_check, $where_to_copy_if_empty) or die('Failed to copy empty VCF: ' . $!);
    } else {
        exec { $commandline_if_nonempty[0] } @commandline_if_nonempty;
        die('Failed to execute command: ' . $!);
    }
}


sub is_compressed {
    my $filename = shift;

    return $filename =~ /\.gz$/;
}

sub vcf_is_empty {
    my $vcf = shift;

    local $SIG{PIPE} = sub { return 1; }; #prevent zgrep "broken pipe" errors

    my $is_compressed = is_compressed($vcf);

    my $cmd = $is_compressed? '/bin/zgrep' : '/bin/grep';

    my $rv = system($cmd, '-q', '-v', '#', $vcf);

    if ($rv == -1) {
        die('Failed to execute command to inspect VCF: ' . $!);
    }

    #0 rv from grep indicates a match for a data line was found => non-empty
    #1 rv from grep indicates only header lines were present => empty
    #other return codes are bad.

    $rv >>= 8;

    if ($rv == 0 or $rv == 1) {
        return $rv;
    } else {
        die('Failed to determine VCF status. Unexpected return code from grep: ' . $rv);
    }
}

sub copy_file {
    my $source = shift;
    my $dest = shift;

    my $source_is_compressed = is_compressed($source);
    my $dest_is_compressed = is_compressed($dest);

    if ($source_is_compressed eq $dest_is_compressed) {
        return copy($source, $dest);
    } elsif ($source_is_compressed) {
        return run(['/bin/gunzip', '-c', $source], '>', $dest);
    } elsif ($dest_is_compressed) {
        return run(['/bin/gzip', '-c', $source], '>', $dest);
    } else {
        die ('Logic error in copy_file.');
    }
}

