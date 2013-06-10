set( ARCHIVE_URL     "https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz" )
set( ARCHIVE_MD5     "0fa47a5ed8fedefab766592785c85ee7" )

set( ARCHIVE_NAME    "yaml-cpp-0.5.1.tar.gz" )
set( ARCHIVE_NAME_WE "yaml-cpp-0.5.1" )

set( PROJECT_CMAKELISTS "
project(YAML_CPP)

set(YAML_CPP_VERSION_MAJOR \"0\")
set(YAML_CPP_VERSION_MINOR \"5\")
set(YAML_CPP_VERSION_PATCH \"1\")
set(YAML_CPP_VERSION \"\${YAML_CPP_VERSION_MAJOR}.\${YAML_CPP_VERSION_MINOR}.\${YAML_CPP_VERSION_PATCH}\")

set(header_directory \"include/yaml-cpp/\")

file(GLOB sources \"src/[a-zA-Z]*.cpp\")
file(GLOB_RECURSE public_headers \"include/yaml-cpp/[a-zA-Z]*.h\")
file(GLOB private_headers \"src/[a-zA-Z]*.h\")
file(GLOB contrib_sources \"src/contrib/[a-zA-Z]*.cpp\")
file(GLOB contrib_public_headers \"include/yaml-cpp/contrib/[a-zA-Z]*.h\")
file(GLOB contrib_private_headers \"src/contrib/[a-zA-Z]*.h\")

include_directories(\${YAML_CPP_SOURCE_DIR}/src)
include_directories(\${YAML_CPP_SOURCE_DIR}/include)

find_package(Boost REQUIRED)
include_directories(\${Boost_INCLUDE_DIRS})

add_library(yaml-cpp STATIC
	\${sources}
	\${public_headers}
	\${private_headers}
	\${contrib_sources}
	\${contrib_public_headers}
	\${contrib_private_headers}
)

set_target_properties(yaml-cpp PROPERTIES
	VERSION \"\${YAML_CPP_VERSION}\"
	SOVERSION \"\${YAML_CPP_VERSION_MAJOR}.\${YAML_CPP_VERSION_MINOR}\"
	PROJECT_LABEL \"yaml-cpp \${LABEL_SUFFIX}\"
)
")