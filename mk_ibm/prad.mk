OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o print_radiation_file.o read_radiation_file.o \
read_radiation_user.o wam_coordinate_module.o

prad:
	mpxlf90 $(OBJECTS) -o prad
