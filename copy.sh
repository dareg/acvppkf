#!/bin/sh

SRC="/home/gmap/mrpm/grassetj/pack/48t2_acvppkf.01.MIMPIIFC1805.2y/acvppkf_extracted/src/"
cp $SRC/acvppkf.F90                                           src/local/arpifs/phys_dmn/acvppkf.F90
cp $SRC/acvppkf.h                                             src/local/arpifs/phys_dmn/acvppkf.h
cp $SRC/convect_chem_transport.F90                            src/local/mpa/conv/internals/convect_chem_transport.F90
cp $SRC/convect_chem_transport.h                              src/local/mpa/conv/internals/convect_chem_transport.h
cp $SRC/convect_closure_adjust_shal.F90                       src/local/mpa/conv/internals/convect_closure_adjust_shal.F90
cp $SRC/convect_closure_adjust_shal.h                         src/local/mpa/conv/internals/convect_closure_adjust_shal.h
cp $SRC/convect_closure_shal.F90                              src/local/mpa/conv/internals/convect_closure_shal.F90
cp $SRC/convect_closure_shal.h                                src/local/mpa/conv/internals/convect_closure_shal.h
cp $SRC/convect_closure_thrvlcl.F90                           src/local/mpa/conv/internals/convect_closure_thrvlcl.F90
cp $SRC/convect_closure_thrvlcl.h                             src/local/mpa/conv/internals/convect_closure_thrvlcl.h
cp $SRC/convect_condens.F90                                   src/local/mpa/conv/internals/convect_condens.F90
cp $SRC/convect_condens.h                                     src/local/mpa/conv/internals/convect_condens.h
cp $SRC/convect_mixing_funct.F90                              src/local/mpa/conv/internals/convect_mixing_funct.F90
cp $SRC/convect_mixing_funct.h                                 src/local/mpa/conv/internals/convect_mixing_funct.h
cp $SRC/convect_satmixratio.F90                               src/local/mpa/conv/internals/convect_satmixratio.F90
cp $SRC/convect_satmixratio.h                                 src/local/mpa/conv/internals/convect_satmixratio.h
cp $SRC/convect_trigger_shal.F90                              src/local/mpa/conv/internals/convect_trigger_shal.F90
cp $SRC/convect_trigger_shal.h                                src/local/mpa/conv/internals/convect_trigger_shal.h
cp $SRC/convect_updraft_shal.F90                              src/local/mpa/conv/internals/convect_updraft_shal.F90
cp $SRC/convect_updraft_shal.h                                src/local/mpa/conv/internals/convect_updraft_shal.h
cp $SRC/ini_convpar.F90                                       src/local/mpa/conv/internals/ini_convpar.h
cp $SRC/modd_convpar.F90                                      src/local/mpa/conv/module/modd_convpar.F90
cp $SRC/modd_convpar_shal.F90                                 src/local/mpa/conv/module/modd_convpar_shal.F90
cp $SRC/modd_convparext.F90                                   src/local/mpa/conv/module/modd_convparext.F90
cp $SRC/modd_cst.F90                                          src/local/mpa/micro/module/modd_cst.F90 
cp $SRC/modd_dimphyexn.F90                                    src/local/mpa/micro/module/modd_dimphyexn.F90 
cp $SRC/modd_nsv.F90                                          src/local/mpa/micro/module/modd_nsv.F90
cp $SRC/shallow_convection.F90                                src/local/mpa/conv/internals/shallow_convection.F90
cp $SRC/shallow_convection.h                                  src/local/mpa/conv/internals/shallow_convection.h
cp $SRC/shallow_convection_compute.F90                        src/local/mpa/conv/internals/shallow_convection_compute.F90
cp $SRC/shallow_convection_compute.h                          src/local/mpa/conv/internals/shallow_convection_compute.h
cp $SRC/shallow_convection_select.F90                         src/local/mpa/conv/internals/shallow_convection_select.F90
cp $SRC/shallow_convection_select.h                           src/local/mpa/conv/internals/shallow_convection_select.h
cp $SRC/apl_arpege_shallow_convection_and_turbulence.F90      src/local/arpifs/phys_dmn/apl_arpege_shallow_convection_and_turbulence.F90
cp $SRC/aplpar.F90                                            src/local/arpifs/phys_dmn/aplpar.F90
cp $SRC/ini_cst.F90                                           src/local/mpa/micro/internals/ini_cst.F90

#For debug only
cp $SRC/yomdbg.F90                                            src/local/mpa/conv/internals/yomdbg.F90
