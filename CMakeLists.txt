cmake_minimum_required(VERSION 2.8)
project( projectZ )
#find_package( vigra REQUIRED )
INCLUDE_DIRECTORIES(/usr/local/include/vigra)

add_executable( projectZ projectZ.cpp )

# target_link_libraries( projectZ ${vigra_LIBS} )
TARGET_LINK_LIBRARIES( projectZ /usr/local/lib/libvigraimpex.so)
