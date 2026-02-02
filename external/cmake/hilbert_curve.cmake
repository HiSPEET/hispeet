add_library(HilbertCurve 
            external/hilbert_curve/hilbert_curve.f 
            )

set_property(TARGET HilbertCurve
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules
             )
