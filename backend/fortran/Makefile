# Makefile for NumBAT. Requires some adjustment to reflect your installation.

COMPILER_VENDOR = gnu95
#COMPILER_VENDOR = intelem

SUBROUTINES_FOR_PYTHON = conv_gmsh.f 

NumBAT.so: NumBAT.pyf *.f
	f2py -c NumBAT.pyf *.f --fcompiler=$(COMPILER_VENDOR) \
	--link-lapack_opt --link-blas_opt

NumBAT.pyf: $(SUBROUTINES_FOR_PYTHON)
	f2py -m NumBAT -h NumBAT.pyf $(SUBROUTINES_FOR_PYTHON) --overwrite-signature

purge:
	-rm NumBAT.pyf
	-rm NumBAT.so
	-rm lib/umf4_f77zwrapper64.o
clean:
	-rm NumBAT.so
