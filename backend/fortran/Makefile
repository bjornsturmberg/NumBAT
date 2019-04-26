# Makefile for NumBAT. May require some adjustment to reflect your installation.

COMPILER_VENDOR = gnu95
#COMPILER_VENDOR = intelem

SUBROUTINES_FOR_PYTHON = conv_gmsh.f py_calc_modes.f \
	py_calc_modes_AC.f EM_mode_energy_int_v2_Ez.f \
	EM_mode_energy_int_Ez.f photoelastic_int.f photoelastic_int_v2.f \
	AC_mode_power_int.f AC_mode_power_int_v2.f AC_mode_power_int_v4.f \
	AC_alpha_int.f AC_alpha_int_v2.f moving_boundary.f array_material_AC.f \
	AC_mode_elastic_energy_int.f AC_mode_elastic_energy_int_v4.f EM_mode_E_energy_int.f \
	H_mode_field_Ez.f array_size.f

LIB_LOCATION = lib
UMFPACK_NAME = umf4_f77zwrapper64.o
#SUITESPARSE_VERSION=3.7.1
SUITESPARSE_VERSION=4.2.0

NumBAT.so: NumBAT.pyf *.f $(LIB_LOCATION)/$(UMFPACK_NAME)
	f2py3 -c NumBAT.pyf *.f --fcompiler=$(COMPILER_VENDOR) \
	--link-lapack_opt --link-blas_opt \
	$(LIB_LOCATION)/$(UMFPACK_NAME) --link-umfpack

lib/umf4_f77zwrapper64.o: lib/umf4_f77zwrapper_SS_$(SUITESPARSE_VERSION).c
	gcc -O5 -DZLONG -I/usr/include/suitesparse -c lib/umf4_f77zwrapper_SS_$(SUITESPARSE_VERSION).c \
		-fPIC -o lib/umf4_f77zwrapper64.o

NumBAT.pyf: $(SUBROUTINES_FOR_PYTHON)
	f2py3 -m NumBAT -h NumBAT.pyf $(SUBROUTINES_FOR_PYTHON) --overwrite-signature

purge:
	-rm NumBAT.so NumBAT.*.so
	-rm lib/umf4_f77zwrapper64.o
clean:
	-rm NumBAT.pyf
	-rm NumBAT.so NumBAT.*.so
