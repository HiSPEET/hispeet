add_library(KdTree 
            external/kdtree/kd_tree.f 
            external/kdtree/kd_tree__precision.f 
            external/kdtree/kd_tree__priority_queue.f 
            )

set_property(TARGET KdTree
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules
             )

target_include_directories(KdTree
                           PUBLIC ${PROJECT_BINARY_DIR}/modules
                           )
