OBJECTS = wam_file_module.o wam_general_module.o wam_timopt_module.o wam_fre_dir_module.o \
wam_interface_module.o wam_grid_module.o wam_current_module.o wam_model_module.o \
wam_ice_module.o wam_output_module.o wam_wind_module.o wam_boundary_module.o \
wam_source_module.o wam_propagation_module.o preproc_module.o wam_coldstart_module.o \
wam_restart_module.o wam_initial_module.o wam_mpi_module.o wam_output_set_up_module.o \
wam_topo_module.o wam_radiation_module.o wam_nest_module.o wam_user_module.o \
wam_special_module.o read_topo_input.o chief.o wavemdl.o initmdl.o read_wam_user.o \
print_wam_status.o read_wind_input.o read_current_input.o wamodel.o read_boundary_input.o \
read_ice_input.o jafu.o wam_mpi_comp_module.o wam_assi_set_up_module.o \
wam_assi_module.o wam_coordinate_module.o readsat.o

chief:
	mpxlf90 $(OBJECTS) -o chief
