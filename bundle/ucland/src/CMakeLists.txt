list ( APPEND _files
src/LND_domains.F90
src/LND_coms.F90
src/LND_io.F90
src/LND_get_input.F90
src/LND_error_handler.F90
src/LND_file_parser.F90
src/LND_string_functions.F90
src/LND_time_manager.F90
src/LND_document.F90
src/version_variable.h
)

set ( src_files
      ${_files}
      PARENT_SCOPE
)
