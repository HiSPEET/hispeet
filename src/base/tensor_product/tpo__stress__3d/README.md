To compare different versions of divergence-penalization, the stress TPO comes in two variants:

1. The default variant, using a physically motivated bulk viscosity in

         tpo__stress__3d_dlci__gen_kernel.f
         tpo__stress__3d_dlci__xsmm_kernel.f

2. The penalty variant, introducing ad-hoc grad-div and flux continuity terms in

         tpo__stress__3d_dlci__gen_kernel-ins_penalty.f
         tpo__stress__3d_dlci__xsmm_kernel-ins_penalty.f

To use the penalty variant, `cmake` must be executed with the option `-D INS_PENALTY=1`. 