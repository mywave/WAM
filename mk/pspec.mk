OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o print_spectra_file.o read_spectra_file.o \
read_spectra_user.o wam_coordinate_module.o wam_output_parameter_module.o 

pspec:
	mpiifort $(OBJECTS) -o pspec
