MACRO( GET_PROJECT PROJ )
  include( ${PROJ} )

  MESSAGE( STATUS "Downloading or verifying ${ARCHIVE_NAME_WE}" )
  FILE( DOWNLOAD ${ARCHIVE_URL}
    "${CMAKE_CURRENT_BINARY_DIR}/${ARCHIVE_NAME}"
    SHOW_PROGRESS STATUS status LOG log EXPECTED_MD5 ${ARCHIVE_MD5}
    )

  list(GET status 0 status_code)
  list(GET status 1 status_string)

  if(NOT status_code EQUAL 0)
    message(FATAL_ERROR "error: downloading '${ARCHIVE_URL}' failed
      status_code: ${status_code}
      status_string: ${status_string}
      log: ${log}
      ")
  endif()

  if( "${CMAKE_CURRENT_BINARY_DIR}/${ARCHIVE_NAME}" IS_NEWER_THAN "${CMAKE_CURRENT_BINARY_DIR}${ARCHIVE_NAME_WE}" )
    MESSAGE( STATUS "Extracting ${ARCHIVE_NAME_WE}" )

    EXECUTE_PROCESS(
      COMMAND "${CMAKE_COMMAND}" -E
      tar -xf "${CMAKE_CURRENT_BINARY_DIR}/${ARCHIVE_NAME}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )

    FILE( WRITE
      "${CMAKE_CURRENT_BINARY_DIR}/${ARCHIVE_NAME_WE}/CMakeLists.txt"
      ${PROJECT_CMAKELISTS}
      )

    ADD_SUBDIRECTORY(
      "${CMAKE_CURRENT_BINARY_DIR}/${ARCHIVE_NAME_WE}"
      )
  endif ()

  unset( ARCHIVE_URL )
  unset( ARCHIVE_MD5 )
  unset( ARCHIVE_NAME )
  unset( ARCHIVE_NAME_WE )
  unset( PROJECT_CMAKELISTS )
ENDMACRO( GET_PROJECT )