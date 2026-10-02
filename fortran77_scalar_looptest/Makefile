FC = gfortran
FFLAGS = -mcmodel=medium -ffixed-line-length-none

exec: torelons_4D_adjoint.f luesher.f
	$(FC) $(FFLAGS) torelons_4D_adjoint.f luesher.f -o exec

clean:
	rm -f exec old_loop_outputs.dat
