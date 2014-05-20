OBJECTS = wam_mpi_module.o wam_file_module.o wam_general_module.o wam_timopt_module.o \
wam_model_module.o wam_source_module.o wam_fre_dir_module.o wam_interface_module.o \
wam_grid_module.o wam_current_module.o wam_special_module.o wam_nest_module.o \
wam_ice_module.o wam_output_module.o wam_print_module.o wam_general_module.o \
wam_output_set_up_module.o wam_mpi_comp_module.o wam_netcdf_module.o wam_coordinate_module.o \
read_current_input.o read_ice_input.o wam_topo_module.o read_topo_input.o jafu.o \
make_netcdf.o 

pnetcdf:
	mpxlf90 $(OBJECTS) -o pnetcdf -qextname -q64 -L/sw/aix61/netcdf-4.1.1-rc1/lib -lnetcdf \
                           -L/sw/aix61/hdf5-1.8.4-patch1/lib -lhdf5 -lhdf5_hl  \
                           -L/sw/aix53/szip-2.1/lib -lsz -L/sw/aix53/zlib-1.2.3/lib -lz -lm
