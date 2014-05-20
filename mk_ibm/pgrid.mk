OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_coordinate_module.o \
wam_print_user_module.o print_grid_file.o read_grid_file.o read_grid_user.o

pgrid:
	mpxlf90 $(OBJECTS) -o pgrid
