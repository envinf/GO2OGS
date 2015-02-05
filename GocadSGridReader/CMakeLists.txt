
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/Utils/FileConverter
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/OpenSG
)

## Gocad SGrid to VTK file converter ##
ADD_EXECUTABLE (GocadSGridReader GocadSGridReaderMain.cpp GocadSGridReader.cpp GocadSGridReader.h GocadNode.h GocadNode.cpp)
SET_TARGET_PROPERTIES(GocadSGridReader PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES(GocadSGridReader
	GeoLib
	MeshLib
	FileIO
	${Boost_LIBRARIES}
)


