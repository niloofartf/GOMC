add_executable(GOMC ${singleExecutable} ${headers} ${libHeaders})
set_target_properties(GOMC PROPERTIES 
  OUTPUT_NAME ${GOMC_name})