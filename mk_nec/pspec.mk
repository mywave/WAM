OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o print_spectra_file.o read_spectra_file.o \
wam_coordinate_module.o read_spectra_user.o

pspec:
	sxmpif90 $(OBJECTS) -o pspec

