project:          HiSPEET
summary:          HiSPEET - High Performance SPEctral Element Techniques.
author:           Jörg Stiller
css:              hispeet.css
src_dir:          ../src/base
exclude_dir:      ../src/base/tensor_product
exclude:          fftw_binding.f
exclude:          boundary_variable__3d-orig.f
exclude:          boundary_variable__3d-alloc.f
exclude:          boundary_variable__3d-pointer.f
exclude:          mesh__3d__mp_alignfrommeshface-body.f
exclude:          mesh__3d__mp_alignwithmeshface-body.f
exclude:          mesh__3d__mp_alignfromneighborface-body.f
exclude:          mesh__3d__mp_alignwithneighborface-body.f
output_dir:       ./html
fixed_extensions: for
                  FOR
extensions:       f
fpp_extensions:   F
preprocess:       false
display:          public
                  protected
                  private
source:           false
graph:            true
graph_maxdepth:   10
graph_nodes:      100
search:           true
docmark:          <
<!--
extra_filetypes: c   //
                 sh  #
                 py  #
-->

