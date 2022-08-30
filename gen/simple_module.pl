#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w

use XML::LibXML;
use strict;

sub gen{
    my $module_name = $_[0];
    my @module_vars = @{$_[1]};

    my $header .= << "EOF";
    MODULE UTIL_$module_name
    USE $module_name
    INTEGER :: CHECK_REF = 1234567890
    CONTAINS
EOF

    my $write_stmts = "";
    foreach my $var (@module_vars){
	$write_stmts .= "\n    WRITE(FH) $var->[1]";
    }

    my $sub_save = "
    SUBROUTINE SAVE_$module_name(FH)
    INTEGER, INTENT(IN) :: FH
    WRITE(FH) CHECK_REF
$write_stmts
    WRITE(FH) CHECK_REF
    END SUBROUTINE SAVE_$module_name\n";

    my $read_stmts = "";
    foreach my $var (@module_vars){
	$read_stmts .= "\n    READ(FH) $var->[1]";
    }

    my $sub_load = "
    SUBROUTINE LOAD_$module_name(FH)
    INTEGER, INTENT(IN) :: FH
    INTEGER :: CHECK
    READ(FH) CHECK
    IF(CHECK /= CHECK_VAL)WRITE(*,*)__LINE__,\"WRONG CONTROL NUMBER\",CHECK
$read_stmts
    READ(FH) CHECK
    IF(CHECK /= CHECK_VAL)WRITE(*,*)__LINE__,\"WRONG CONTROL NUMBER\",CHECK
    END SUBROUTINE LOAD_$module_name\n";

    my $footer = "END MODULE UTIL_$module_name";

    print($header);
    print($sub_save);
    print($sub_load);
    print($footer);
}

my @module_vars = ();
my $module_name = "noname";
my $xml_file = $ARGV[0];

my $uri = 'http://fxtran.net/#syntax';
my $xpc = 'XML::LibXML::XPathContext'->new ();
$xpc->registerNs (fx => $uri);
my $doc = 'XML::LibXML'->load_xml (location => $xml_file);

foreach my $stmt ($xpc->findnodes('.//fx:T-decl-stmt', $doc)){
    my $var_type = $xpc->findnodes('.//fx:T-N', $stmt);
    foreach my $var ($xpc->findnodes('.//fx:EN-decl', $stmt)){
        my $var_name = $xpc->findnodes('.//fx:n', $var);
        push(@module_vars, [$var_type, $var_name]);
    }
}

($module_name) = $xpc->findnodes('//fx:module-stmt/fx:module-N/fx:N/fx:n', $doc);
$module_name = $module_name->to_literal();

gen($module_name, \@module_vars)
