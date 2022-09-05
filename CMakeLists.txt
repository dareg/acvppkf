cmake_minimum_required(VERSION 3.15)


project(acvppkf_extract)
enable_language (Fortran)
enable_language (C)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)


find_package(OpenMP)
if(NOT OpenMP_Fortran_FOUND)
	MESSAGE(FATAL_ERROR "Could not find OpenMP support.")
endif()

#add_compile_options(-ffree-form -ffree-line-length-512 -g -fdefault-real-8 -fconvert=big-endian -fcheck=all -O3)

set(SRC
	src/acvppkf.F90
	src/convect_chem_transport.F90
	src/convect_closure_adjust_shal.F90
	src/convect_closure_shal.F90
	src/convect_closure_thrvlcl.F90
	src/convect_condens.F90
	src/convect_mixing_funct.F90
	src/convect_satmixratio.F90
	src/convect_trigger_shal.F90
	src/convect_updraft_shal.F90
	src/convection_shal.F90
	src/crc.c
	src/eint_mod.F90
	src/ini_convpar.F90
	src/lfa_R8I4.F90
	src/modd_convpar.F90
	src/modd_convpar_shal.F90
	src/modd_convparext.F90
	src/modd_cst.F90
	src/modd_nsv.F90
	src/modd_parameters.F90
	src/model_physics_mf_mod.F90
	src/parkind1.F90
	src/posnam.F90
	src/posname.F90
	src/shallow_convection.F90
	src/util_modd_convpar_shal.F90
	src/util_modd_cst.F90
	src/util_model_physics_mf_type_mod.F90
	src/util_sl_struct_mod.F90
	src/util_tarphy_mod.F90
	src/util_tcst_mod.F90
	src/util_tcvmnh_mod.F90
	src/util_tdprecips_mod.F90
	src/util_tdvisi_mod.F90
	src/util_tlouis_mod.F90
	src/util_tmse_mod.F90
	src/util_tnorgwd_mod.F90
	src/util_tparar_mod.F90
	src/util_tphy0_mod.F90
	src/util_tphy1_mod.F90
	src/util_tphy2_mod.F90
	src/util_tphy3_mod.F90
	src/util_tphy_mod.F90
	src/util_tphyds_mod.F90
	src/util_tsimphl_mod.F90
	src/util_ttoph_mod.F90
	src/util_tvdoz_mod.F90
	src/wrap_acvppkf.F90
	src/wrscmr.F90
	src/yomarphy.F90
	src/yomcst.F90
	src/yomct3.F90
	src/yomcvmnh.F90
	src/yomdbg.F90
	src/yomdprecips.F90
	src/yomdvisi.F90
	src/yomdyncore.F90
	src/yomhook.F90
	src/yomlouis.F90
	src/yomlsforc.F90
	src/yomlun.F90
	src/yomlun_ifsaux.F90
	src/yommp0.F90
	src/yommse.F90
	src/yomnorgwd.F90
	src/yomparar.F90
	src/yomphy.F90
	src/yomphy0.F90
	src/yomphy1.F90
	src/yomphy2.F90
	src/yomphy3.F90
	src/yomphyds.F90
	src/yomsimphl.F90
	src/yomtoph.F90
	src/yomvdoz.F90
	)

add_executable(main.x main.F90 ${SRC})
target_link_libraries(main.x PUBLIC OpenMP::OpenMP_Fortran)
target_compile_options(main.x PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-convert big_endian -assume byterecl -align array64byte,all -traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model source -qopt-report=5 -qopt-report-phase=vec -ftz -diag-disable=remark,cpu-dispatch -free  -DLINUX -DLITTLE_ENDIAN -DLITTLE -DADDRESS64 -DGRIB_API_1 -g -O2 -march=core-avx2 -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma -list -show nomap -r8 -heap-arrays>)
