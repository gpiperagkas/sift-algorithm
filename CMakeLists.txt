cmake_minimum_required(VERSION 2.8)

project( colordescr )

find_package( OpenCV REQUIRED )

include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable( colordescriptor main.cpp  imagergb.cpp featurevector.cpp classifier.cpp)

target_link_libraries( colordescriptor ${OpenCV_LIBS} )