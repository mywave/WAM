OBJECTS = wam_mpi_module.o wam_file_module.o wam_boundary_module.o wam_coordinate_module.o \
wam_general_module.o wam_timopt_module.o wam_fre_dir_module.o \
wam_interface_module.o wam_grid_module.o preproc_module.o \
preproc_user_module.o wam_nest_module.o read_boundary_input.o \
wam_special_module.o wam_output_set_up_module.o \
preproc.o read_topography.o read_preproc_user.o wam_mpi_comp_module.o

preproc:
	sxmpif90 $(OBJECTS) -o preproc

