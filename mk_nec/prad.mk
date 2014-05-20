OBJECTS =  wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o print_radiation_file.o read_radiation_file.o \
wam_coordinate_module.o read_radiation_user.o

prad:
	sxmpif90 $(OBJECTS) -o prad
