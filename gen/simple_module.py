#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import sys

def gen(module_name, module_vars):
    header = \
    f'''    MODULE UTIL_{module_name}
    USE {module_name}
    INTEGER :: CHECK_REF = 1234567890
    CONTAINS'''

    write_stmts = ""
    for var_type, var_name in module_vars:
        write_stmts += f'''\n    WRITE(FH) {var_name}'''

    sub_save = f'''
    SUBROUTINE SAVE_{module_name}(FH)
    INTEGER, INTENT(IN) :: FH
    WRITE(FH) CHECK_REF
{write_stmts}
    WRITE(FH) CHECK_REF
    END SUBROUTINE SAVE_{module_name}'''

    read_stmts = ""
    for var_type, var_name in module_vars:
        read_stmts += f'''\n    READ(FH) {var_name}'''
    
    sub_load = f'''
    SUBROUTINE LOAD_{module_name}(FH)
    INTEGER, INTENT(IN) :: FH
    INTEGER :: CHECK
    READ(FH) CHECK
    IF(CHECK /= CHECK_VAL)WRITE(*,*)__LINE__,\"WRONG CONTROL NUMBER\",CHECK
{read_stmts}
    READ(FH) CHECK
    IF(CHECK /= CHECK_VAL)WRITE(*,*)__LINE__,\"WRONG CONTROL NUMBER\",CHECK
    END SUBROUTINE LOAD_{module_name}'''

    footer=f'''END MODULE UTIL_{module_name}'''

    print(header)
    print(sub_save)
    print(sub_load)
    print(footer)

module_vars = []
module_name="noname"
xml_file = sys.argv[1]

namespaces = {'fx': 'http://fxtran.net/#syntax'}
tree = ET.parse(xml_file)
root = tree.getroot()

for stmt in root.findall('.//fx:T-decl-stmt', namespaces):
    var_type = stmt.find('.//fx:T-N', namespaces).text
    for var in stmt.findall('.//fx:EN-decl', namespaces):
        var_name = var.find('.//fx:n',namespaces).text
        module_vars.append((var_type, var_name))

module_stmt = root.find('.//fx:module-stmt', namespaces)
module_name=module_stmt.find('.//fx:n',namespaces).text

gen(module_name,module_vars)

