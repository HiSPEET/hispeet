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

//> summary:  Wrappers for accessing ParMETIS from Fortran
//> author:   Joerg Stiller
//> date:     2014/10/08
//>
//> The routines provided with this file are based on the orignal wrappers
//> with ParMETIS 4. They are supplemented with distributed explicit Fortran
//> interfaces defined in parmetis_binding.f. Relying on standardized ISO C
//> binding these interfaces provide a safe access to the corresponding ParMETIS
//> routines.
//==============================================================================

#include <mpi.h>
#include <parmetis.h>

//------------------------------------------------------------------------------
//> ParMETIS_V3_PartKway Fortran to C wrapper

void ParMETIS_V3_PartKway_F2C(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
    idx_t *vwgt, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon,
    idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options,
    idx_t *edgecut, idx_t *part, MPI_Fint *comm)
{
  MPI_Comm ccomm = MPI_Comm_f2c(*comm);

  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag,
      ncon, nparts, tpwgts, ubvec, options, edgecut, part, &ccomm);

  return;
}
