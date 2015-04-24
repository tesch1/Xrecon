cmake_minimum_required (VERSION 2.8)

project (Xrecon C)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -std=gnu99")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions (-DDEBUG)
endif ()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set (USE_FFTWD 1)
find_package (FFTW3 REQUIRED)
find_package (GSL REQUIRED)
#find_package (BLAS)
find_package (TIFF REQUIRED)

include_directories (${FFTW3_INCLUDE_DIR})
include_directories (${GSL_INCLUDE_DIR})
include_directories (${TIFF_INCLUDE_DIR})

set (Xrecon_SRCS
  Xrecon.c
  recon1D/recon1D.c
  recon1D/default1D.c
  recon2D/recon2D.c
  recon2D/default2D.c
  recon2D/profile2D.c
  recon3D/recon3D.c
  recon3D/default3D.c
  reconEPI/reconEPI.c
  reconEPI/defaultEPI.c
  reconEPI/dprocEPI.c
  reconEPI/prescanEPI.c
  nifti/niftiwrite.c
  common1D/dread1D.c
  common1D/dproc1D.c
  common1D/dutils1D.c
  common1D/write1D.c
  common2D/dread2D.c
  common2D/dproc2D.c
  common2D/noise2D.c
  common2D/dmask2D.c
  common2D/dutils2D.c
  common2D/fdfwrite2D.c
  common2D/rawwrite2D.c
  common2D/tifwrite2D.c
  common3D/fdfwrite3D.c
  common3D/rawIO3D.c
  common3D/dproc3D.c
  common/dproc.c
  common/dhead.c
  common/dutils.c
  common/options.c
  common/pars.c
  common/utils.c)

add_executable (Xrecon ${Xrecon_SRCS})
target_link_libraries (Xrecon m)
target_link_libraries (Xrecon ${FFTW3_LIBRARY})
target_link_libraries (Xrecon ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
#target_link_libraries (Xrecon ${BLAS_LIBRARIES})
target_link_libraries (Xrecon ${TIFF_LIBRARY})
