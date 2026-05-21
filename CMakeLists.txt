# CMakeList.txt : CMake project for BWRSflow, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("BWRSflow")

# Add source to this project's executable.
add_executable (BWRSflow
    "App/BWRSflow.cpp" "App/BWRSflow.h" "App/BWRSflow.rc"
    "Src/Common/constants.h"
    "Src/Common/types.h"
    "Src/Data/components.h"
    "Src/Data/components.cpp"
    "Src/Formulas/bwrs.h"
    "Src/Formulas/bwrs.cpp"
    "Src/Formulas/viscosity.h"
    "Src/Formulas/viscosity.cpp"
    "Src/Standards/iso5167.h"
    "Src/Standards/iso5167.cpp"
    "Src/Ui/console.h"
    "Src/Ui/console.cpp"
    "Src/Tests/validation.h"
    "Src/Tests/validation.cpp"
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET BWRSflow PROPERTY CXX_STANDARD 20)
endif()

if (MSVC)
  target_compile_definitions(BWRSflow PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# TODO: Add tests and install targets if needed.
