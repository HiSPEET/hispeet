file(COPY cylinder_2d.msh DESTINATION .)
file(COPY pipe.msh        DESTINATION .)

file(ARCHIVE_EXTRACT 
         INPUT ${CMAKE_SOURCE_DIR}/develop/gmsh/linear_cascade.msh.zip 
         DESTINATION .)
