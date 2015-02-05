INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(resample resample.cpp)
TARGET_LINK_LIBRARIES(resample ${VTK_LIBRARIES})

ADD_EXECUTABLE(IntegrateFaultFacesets IntegrateFaultFacesets.cpp)
TARGET_LINK_LIBRARIES(IntegrateFaultFacesets ${VTK_LIBRARIES})
