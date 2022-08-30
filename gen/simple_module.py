def gen(module_name, module_vars):
    header = f'''
    MODULE UTIL_{module_name}
    USE {module_name}
    CONTAINS'''

    write_stmts = ""
    for var_type, var_name in module_vars:
        write_stmts += f'''
    WRITE(FH) {var_name}'''

    sub_save = f'''
    SUBROUTINE SAVE_{module_name}(FH)
    INTEGER, INTENT(IN) :: FH
    {write_stmts}
    END SUBROUTINE SAVE_{module_name}
    '''

    read_stmts = ""
    for var_type, var_name in module_vars:
        read_stmts += f'''
    READ(FH) {var_name}'''
    
    sub_load = f'''
    SUBROUTINE LOAD_{module_name}(FH)
    INTEGER, INTENT(IN) :: FH
    {read_stmts}
    END SUBROUTINE LOAD_{module_name}
    '''

    footer=f'''
    END MODULE UTIL_{module_name}
    '''

    print(header)
    print(sub_save)
    print(sub_load)
    print(footer)

import xml.etree.ElementTree as ET
import sys

module_vars = []
module_name="noname"
xml_file = sys.argv[1]

namespaces = {'fx': 'http://fxtran.net/#syntax'}
tree = ET.parse(xml_file)
root = tree.getroot()

for stmt in root.findall('.//fx:T-decl-stmt', namespaces):
    var_type = stmt.find('.//fx:T-N', namespaces).text
    var_name = stmt.find('.//fx:n',namespaces).text
    module_vars.append((var_type, var_name))

module_stmt = root.find('.//fx:module-stmt', namespaces)
module_name=module_stmt.find('.//fx:n',namespaces).text

gen(module_name,module_vars)

