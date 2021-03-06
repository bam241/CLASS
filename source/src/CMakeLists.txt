message( STATUS "Building Dict files from root Cint... ")

ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/CLASSObjectDict ${CLASS_HEADERS_DIR_INCLUDE}/CLASSObject.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/CLASSFacilityDict ${CLASS_HEADERS_DIR_INCLUDE}/CLASSFacility.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/CLASSBackEndDict ${CLASS_HEADERS_DIR_INCLUDE}/CLASSBackEnd.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/StorageDict ${CLASS_HEADERS_DIR_INCLUDE}/Storage.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/ReactorDict ${CLASS_HEADERS_DIR_INCLUDE}/Reactor.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/FabricationPlantDict ${CLASS_HEADERS_DIR_INCLUDE}/FabricationPlant.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/PoolDict ${CLASS_HEADERS_DIR_INCLUDE}/Pool.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/SeparationPlantDict ${CLASS_HEADERS_DIR_INCLUDE}/SeparationPlant.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/IsotopicVectorDict ${CLASS_HEADERS_DIR_INCLUDE}/IsotopicVector.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/ZAIDict ${CLASS_HEADERS_DIR_INCLUDE}/ZAI.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/EvolutionDataDict ${CLASS_HEADERS_DIR_INCLUDE}/EvolutionData.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/DecayDataBankDict ${CLASS_HEADERS_DIR_INCLUDE}/DecayDataBank.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/ScenarioDict ${CLASS_HEADERS_DIR_INCLUDE}/Scenario.hxx )
ROOT_GENERATE_DICTIONARY( ${CMAKE_CURRENT_SOURCE_DIR}/ReactorSchedulerDict ${CLASS_HEADERS_DIR_INCLUDE}/ReactorScheduler.hxx )

add_library(CLASSpkg_root STATIC CLASSBackEnd.cxx CLASSBackEndDict.cxx CLASSObject.cxx CLASSObjectDict.cxx CLASSFacility.cxx CLASSFacilityDict.cxx Storage.cxx StorageDict.cxx Reactor.cxx ReactorDict.cxx FabricationPlant.cxx FabricationPlantDict.cxx Pool.cxx PoolDict.cxx SeparationPlant.cxx SeparationPlantDict.cxx IsotopicVector.cxx IsotopicVectorDict.cxx ZAI.cxx ZAIDict.cxx EvolutionData.cxx EvolutionDataDict.cxx DecayDataBank.cxx DecayDataBankDict.cxx Scenario.cxx ScenarioDict.cxx ReactorScheduler.cxx ReactorSchedulerDict.cxx)

target_link_libraries(CLASSpkg_root ${ROOT_LIBRARIES})

message("------- CLASSpkg_root library DONE... -------")

LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSLogger.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ZAI.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/IsotopicVector.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ZAIMass.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ZAIHeat.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ZAITox.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSNucleiFiliation.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSObject.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSFacility.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/FabricationPlant.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Reactor.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSBackEnd.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/SeparationPlant.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Storage.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Pool.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/DecayDataBank.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/DynamicalSystem.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/IrradiationModel.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/EquivalenceModel.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/XSModel.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/PhysicsModels.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/EvolutionData.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReactorScheduler.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Scenario.cxx")
LIST( APPEND CLASS_SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CLASSReader.cxx")

SET( CLASS_SRC_FILE ${CLASS_SRC_FILE} PARENT_SCOPE)

message("------- CLASS_SRC_FILE VARIABLE -------")
MESSAGE( STATUS "CLASS_SRC_FILE=" ${CLASS_SRC_FILE})
