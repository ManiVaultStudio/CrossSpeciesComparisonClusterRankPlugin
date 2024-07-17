cmake_minimum_required(VERSION 3.17)

# -----------------------------------------------------------------------------
# CrossSpeciesComparisonPhyTVu Plugin
# -----------------------------------------------------------------------------
set(CrossSpeciesComparisonPhyTVu "CrossSpeciesComparisonClusterRankPlugin")
PROJECT(${CrossSpeciesComparisonPhyTVu})

# -----------------------------------------------------------------------------
# CMake Options
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DWIN32 /EHsc /MP /permissive- /Zc:__cplusplus")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif(MSVC)

# -----------------------------------------------------------------------------
# Set install directory
# -----------------------------------------------------------------------------
# Check if the directory to the ManiVault installation has been provided
if(NOT DEFINED MV_INSTALL_DIR)
    set(MV_INSTALL_DIR "" CACHE PATH "Directory where ManiVault is installed")
    message(FATAL_ERROR "Please set MV_INSTALL_DIR to the directory where ManiVault is installed")
endif()
file(TO_CMAKE_PATH ${MV_INSTALL_DIR} MV_INSTALL_DIR)

if(NOT DEFINED MV_CSCTD_INSTALL_DIR) 
    set(MV_CSCTD_INSTALL_DIR ${MV_INSTALL_DIR} CACHE PATH "Directory where TreeData is installed")
endif()

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Widgets WebEngineWidgets Concurrent REQUIRED)

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
# Define the plugin sources
set(CrossSpeciesComparisonPhyTVu_SOURCES
    src/CrossSpeciesComparisonClusterRankPlugin.h
    src/CrossSpeciesComparisonClusterRankPlugin.cpp
    src/ChartWidget.h
    src/ChartWidget.cpp
	src/ChartOptions.h
    src/ChartOptions.cpp
    src/CrossSpeciesComparisonClusterRankPlugin.json
)

set(PLUGIN_MOC_HEADERS
	src/CrossSpeciesComparisonClusterRankPlugin.h
)

set(WEB_MAIN
    res/CrossSpeciesComparisonPhyTVu_chart/CrossSpeciesComparisonPhyTVu_chart.html
    res/CrossSpeciesComparisonPhyTVu_chart/CrossSpeciesComparisonPhyTVu_chart.js
    res/CrossSpeciesComparisonPhyTVu_chart/ChartFile.js
    res/CrossSpeciesComparisonPhyTVu_chart/DataFile.js
    res/CrossSpeciesComparisonPhyTVu_chart/FunctionFile.js
	res/CrossSpeciesComparisonPhyTVu_chart/CrossSpeciesComparisonPhyTVu_chart.css
)

set(WEB_LIBS
    res/js_libs/d3.v6.js
    res/js_libs/d3.v4.js
    res/js_libs/d3.v5.js
    res/js_libs/d3-scale-chromatic.v1.min.js
    res/js_libs/qwebchannel.tools.js
)

set(Lib
    src/lib/Clustering/fastcluster.cpp
	src/lib/Clustering/fastcluster.h
	src/lib/Clustering/fastcluster_dm.h
	src/lib/Clustering/fastcluster_R_dm.h
	src/lib/JSONnlohmann/json.hpp
	src/lib/JSONnlohmann/json_fwd.hpp
)
set(AUX
    res/CrossSpeciesComparisonPhyTVu_chart.qrc
)

qt6_add_resources(RESOURCE_FILES res/CrossSpeciesComparisonPhyTVu_chart.qrc)

source_group( Plugin FILES ${CrossSpeciesComparisonPhyTVu_SOURCES})
source_group(Web\\MainCode FILES ${WEB_MAIN})
source_group(Web\\LibraryFiles FILES ${WEB_LIBS})
source_group( ResourceFiles FILES ${RESOURCES})
source_group( LibFiles FILES ${Lib})
source_group( ResourceFiles\\generated FILES ${CMAKE_CURRENT_BINARY_DIR}/qrc_phylo_tree_view_resources.cpp)
source_group( Aux FILES ${AUX})
# -----------------------------------------------------------------------------
# CMake Target
# -----------------------------------------------------------------------------
# Create dynamic library for the plugin
add_library(${CrossSpeciesComparisonPhyTVu} SHARED ${CrossSpeciesComparisonPhyTVu_SOURCES} ${AUX} ${RESOURCES} ${RESOURCE_FILES} ${WEB_MAIN} ${WEB_LIBS} ${Lib} )
qt_wrap_cpp(CrossSpeciesComparisonPhyTVu_MOC ${PLUGIN_MOC_HEADERS} TARGET ${CrossSpeciesComparisonPhyTVu})
target_sources(${CrossSpeciesComparisonPhyTVu} PRIVATE ${CrossSpeciesComparisonPhyTVu_MOC})

# -----------------------------------------------------------------------------
# Target include directories
# -----------------------------------------------------------------------------
# Include HDPS headers, including system data plugins
target_include_directories(${CrossSpeciesComparisonPhyTVu} PRIVATE "${MV_INSTALL_DIR}/$<CONFIGURATION>/include/")
target_include_directories(${CrossSpeciesComparisonPhyTVu} PRIVATE "${MV_CSCTD_INSTALL_DIR}/$<CONFIGURATION>/include/")

# -----------------------------------------------------------------------------
# Target properties
# -----------------------------------------------------------------------------
# Request C++17
target_compile_features(${CrossSpeciesComparisonPhyTVu} PRIVATE cxx_std_17)

# -----------------------------------------------------------------------------
# Target library linking
# -----------------------------------------------------------------------------
# Link to Qt libraries
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE Qt6::Widgets)
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE Qt6::WebEngineWidgets)
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE Qt6::Concurrent)

# Link to HDPS and data plugins
# The link path in this repo assumes that the ManiVault core was built locally
# in contrast to having been installed with an installer. In the latter case you'll have
# to adapt the MV_LINK_PATH and PLUGIN_LINK_PATH to your install folder
set(MV_LINK_PATH "${MV_INSTALL_DIR}/$<CONFIGURATION>/lib")
set(PLUGIN_LINK_PATH "${MV_INSTALL_DIR}/$<CONFIGURATION>/$<IF:$<CXX_COMPILER_ID:MSVC>,lib,Plugins>")
set(CSCTDPLUGIN_LINK_PATH "${MV_CSCTD_INSTALL_DIR}/$<CONFIGURATION>/$<IF:$<CXX_COMPILER_ID:MSVC>,lib,Plugins>")
set(MV_LINK_SUFFIX $<IF:$<CXX_COMPILER_ID:MSVC>,${CMAKE_LINK_LIBRARY_SUFFIX},${CMAKE_SHARED_LIBRARY_SUFFIX}>)

set(MV_LINK_LIBRARY "${MV_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}MV_Public${MV_LINK_SUFFIX}")
set(POINTDATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}PointData${MV_LINK_SUFFIX}") 
set(CLUSTERDATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}ClusterData${MV_LINK_SUFFIX}") 
set(CROSSSPECIESCOMPARISONTREEDATA_LINK_LIBRARY "${CSCTDPLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}CrossSpeciesComparisonTreeData${MV_LINK_SUFFIX}") 

target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE "${MV_LINK_LIBRARY}")
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE "${POINTDATA_LINK_LIBRARY}")
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE "${CLUSTERDATA_LINK_LIBRARY}")
target_link_libraries(${CrossSpeciesComparisonPhyTVu} PRIVATE "${CROSSSPECIESCOMPARISONTREEDATA_LINK_LIBRARY}")

# -----------------------------------------------------------------------------
# Target installation
# -----------------------------------------------------------------------------
# Install the shared plugin libary to the "Plugins" folder in the ManiVault install directory
install(TARGETS ${CrossSpeciesComparisonPhyTVu}
    RUNTIME DESTINATION Plugins COMPONENT PLUGINS # Windows .dll
    LIBRARY DESTINATION Plugins COMPONENT PLUGINS # Linux/Mac .so
)

add_custom_command(TARGET ${CrossSpeciesComparisonPhyTVu} POST_BUILD
    COMMAND "${CMAKE_COMMAND}"
        --install ${CMAKE_CURRENT_BINARY_DIR}
        --config $<CONFIGURATION>
        --prefix ${MV_INSTALL_DIR}/$<CONFIGURATION>
)

# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------
# Automatically set the debug environment (command + working directory) for MSVC
if(MSVC)
    set_property(TARGET ${CrossSpeciesComparisonPhyTVu} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY $<IF:$<CONFIG:DEBUG>,${MV_INSTALL_DIR}/debug,${MV_INSTALL_DIR}/release>)
    set_property(TARGET ${CrossSpeciesComparisonPhyTVu} PROPERTY VS_DEBUGGER_COMMAND $<IF:$<CONFIG:DEBUG>,"${MV_INSTALL_DIR}/debug/ManiVault Studio.exe","${MV_INSTALL_DIR}/release/ManiVault Studio.exe">)
endif()
