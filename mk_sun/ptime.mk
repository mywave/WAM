OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o print_time.o read_time_user.o read_grid_file.o \
wam_coordinate_module.o

ptime:
	mpif90 $(OBJECTS) -o ptime
