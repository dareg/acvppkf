#!/bin/sh

SRC="/home/gmap/mrpm/grassetj/pack/48t2_acvppkf.01.MIMPIIFC1805.2y/acvppkf_extracted/src/"
copy_if_different() {
	cmp -s $1 $2 || cp $1 $2
}
copy_if_different $SRC/acvppkf.F90                                           src/local/arpifs/phys_dmn/acvppkf.F90
copy_if_different $SRC/acvppkf.h                                             src/local/arpifs/phys_dmn/acvppkf.h
copy_if_different $SRC/convect_chem_transport.F90                            src/local/mpa/conv/internals/convect_chem_transport.F90
copy_if_different $SRC/convect_chem_transport.h                              src/local/mpa/conv/internals/convect_chem_transport.h
copy_if_different $SRC/convect_closure_adjust_shal.F90                       src/local/mpa/conv/internals/convect_closure_adjust_shal.F90
copy_if_different $SRC/convect_closure_adjust_shal.h                         src/local/mpa/conv/internals/convect_closure_adjust_shal.h
copy_if_different $SRC/convect_closure_shal.F90                              src/local/mpa/conv/internals/convect_closure_shal.F90
copy_if_different $SRC/convect_closure_shal.h                                src/local/mpa/conv/internals/convect_closure_shal.h
copy_if_different $SRC/convect_closure_thrvlcl.F90                           src/local/mpa/conv/internals/convect_closure_thrvlcl.F90
copy_if_different $SRC/convect_closure_thrvlcl.h                             src/local/mpa/conv/internals/convect_closure_thrvlcl.h
copy_if_different $SRC/convect_condens.F90                                   src/local/mpa/conv/internals/convect_condens.F90
copy_if_different $SRC/convect_condens.h                                     src/local/mpa/conv/internals/convect_condens.h
copy_if_different $SRC/convect_mixing_funct.F90                              src/local/mpa/conv/internals/convect_mixing_funct.F90
copy_if_different $SRC/convect_mixing_funct.h                                src/local/mpa/conv/internals/convect_mixing_funct.h
copy_if_different $SRC/convect_satmixratio.F90                               src/local/mpa/conv/internals/convect_satmixratio.F90
copy_if_different $SRC/convect_satmixratio.h                                 src/local/mpa/conv/internals/convect_satmixratio.h
copy_if_different $SRC/convect_trigger_shal.F90                              src/local/mpa/conv/internals/convect_trigger_shal.F90
copy_if_different $SRC/convect_trigger_shal.h                                src/local/mpa/conv/internals/convect_trigger_shal.h
copy_if_different $SRC/convect_updraft_shal.F90                              src/local/mpa/conv/internals/convect_updraft_shal.F90
copy_if_different $SRC/convect_updraft_shal.h                                src/local/mpa/conv/internals/convect_updraft_shal.h
copy_if_different $SRC/ini_convpar.F90                                       src/local/mpa/conv/internals/ini_convpar.h
copy_if_different $SRC/modd_convpar.F90                                      src/local/mpa/conv/module/modd_convpar.F90
copy_if_different $SRC/modd_convpar_shal.F90                                 src/local/mpa/conv/module/modd_convpar_shal.F90
copy_if_different $SRC/modd_convparext.F90                                   src/local/mpa/conv/module/modd_convparext.F90
copy_if_different $SRC/modd_cst.F90                                          src/local/mpa/micro/module/modd_cst.F90 
copy_if_different $SRC/modd_dimphyexn.F90                                    src/local/mpa/micro/module/modd_dimphyexn.F90 
copy_if_different $SRC/modd_nsv.F90                                          src/local/mpa/micro/module/modd_nsv.F90
copy_if_different $SRC/shallow_convection.F90                                src/local/mpa/conv/internals/shallow_convection.F90
copy_if_different $SRC/shallow_convection.h                                  src/local/mpa/conv/internals/shallow_convection.h
copy_if_different $SRC/shallow_convection_compute.F90                        src/local/mpa/conv/internals/shallow_convection_compute.F90
copy_if_different $SRC/shallow_convection_compute.h                          src/local/mpa/conv/internals/shallow_convection_compute.h
copy_if_different $SRC/shallow_convection_select.F90                         src/local/mpa/conv/internals/shallow_convection_select.F90
copy_if_different $SRC/shallow_convection_select.h                           src/local/mpa/conv/internals/shallow_convection_select.h
copy_if_different $SRC/apl_arpege_shallow_convection_and_turbulence.F90      src/local/arpifs/phys_dmn/apl_arpege_shallow_convection_and_turbulence.F90
copy_if_different $SRC/aplpar.F90                                            src/local/arpifs/phys_dmn/aplpar.F90
copy_if_different $SRC/ini_cst.F90                                           src/local/mpa/micro/internals/ini_cst.F90

#For debug only
copy_if_different $SRC/yomdbg.F90                                            src/local/mpa/conv/internals/yomdbg.F90
