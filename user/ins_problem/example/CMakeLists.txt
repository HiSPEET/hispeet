#===============================================================================
# User-defined INS problem: example
#
# Included by src/ins_flow/CMakeLists.txt when HiSPEET is configured with
# -D INS_USER_PROBLEM=example. See ../README.md.

# user libraries, if needed, preferably as imported targets, e.g.
#   find_package(MyLib REQUIRED)

hispeet_ins_user_problem(
    SOURCES   ins__problem__user__3d.f
  # LIBRARIES MyLib::MyLib
)
