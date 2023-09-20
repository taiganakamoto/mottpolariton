cmake_minimum_required(VERSION 3.20)
project(mottpolariton NONE)
find_program(Julia_EXECUTABLE julia)
message(STATUS "Found Julia executable: ${Julia_EXECUTABLE}")
message(STATUS "project name is automaically set as: ${PROJECT_NAME}")
set (execute_name "main")
enable_testing()

configure_file(${execute_name} ${execute_name} FILE_PERMISSIONS OWNER_READ 
OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
WORLD_READ WORLD_EXECUTE)