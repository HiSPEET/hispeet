//> summary:   Wrappers for accessing METIS from Fortran
//> author:    Joerg Stiller
//> date:      2025/06/07
//> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
//>
//>### Wrappers for accessing METIS from Fortran
//>
//> Relying on standardized ISO C binding these interfaces provide a safe access
//> to the corresponding METIS routines.
//==============================================================================

#include <stddef.h>
#include <metis.h>

//------------------------------------------------------------------------------
//> METIS_PartGraphRecursive Fortran to C wrapper

void METIS_PartGraphRecursive_F2C(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
    idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts,
    real_t *ubvec, idx_t *edgecut, idx_t *part)

{
  idx_t options[METIS_NOPTIONS];

  METIS_SetDefaultOptions(options);

  METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, NULL, adjwgt,
      nparts, tpwgts, ubvec, options, edgecut, part);

  return;

}

//------------------------------------------------------------------------------
//> METIS_PartGraphKway Fortran to C wrapper

void METIS_PartGraphKway_F2C(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
    idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts,
    real_t *ubvec, idx_t *edgecut, idx_t *part)

{
  idx_t options[METIS_NOPTIONS];

  METIS_SetDefaultOptions(options);

  // minimize edge cut (and not total communication volume)
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

  METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, NULL, adjwgt, nparts,
      tpwgts, ubvec, options, edgecut, part);

  return;
}
