#
# Define variables of module groups for use by client projects.
#
macro(setup_module_variables_for_cater_helpers caterhelp_directory maths_directory)

  # load the macro 'setup_module_variables_for_maths' from module_definitions.cmake (found in the maths_directory)
  include(${maths_directory}/cmake/module_definitions.cmake)
  # Use the macro to load the SM_*_MODULES variables
  setup_module_variables_for_maths (${maths_directory} "")

  set(CATERHELP_MATHS_MODULES
    ${SM_VEC_MODULES}
    ${SM_VVEC_MODULES}
    ${SM_ALGO_MODULES}
  )
  list(REMOVE_DUPLICATES CATERHELP_MATHS_MODULES)

  set(CATERHELP_MODULES
    ${CATERHELP_MATHS_MODULES}
    ${caterhelp_directory}/modules/cater_helpers.cppm
  )
  list(REMOVE_DUPLICATES CATERHELP_MODULES)

endmacro()
