//------------------------------------------------------------------------------
// This file is part of HiSPEET: High-order Spectral Element Techniques
//
// Copyright (C) 2026 by the HiSPEET authors and the
// Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
//
// HiSPEET is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HiSPEET is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------------------

//> summary:  Wrappers for accessing METIS from Fortran
//> author:   Joerg Stiller
//> date:     2025/06/07
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
