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

//> summary:  Wrappers for accessing the VTK XML writer from Fortran
//> author:   Joerg Stiller
//> date:     2014/05/16
//>
//> The module provides access to the C routines defined in vtkXMLWriterC.h. The
//> corresponding Fortran interfaces are found in vtk_binding.f.
//> This implementation is derived from vtkXMLWriterF.c authored by Ken Martin,
//> Will Schroeder, Bill Lorensen. In contrast to the original version, only one
//> wrapper is needed for each routine. The wrappers are accessed from Fortran
//> using the standardized ISO C binding.

/* Calls will be forwarded to the C interface.  */
#include "vtkXMLWriterC.h"

#include <stdio.h>  /* fprintf  */

/* Define a static-storage default-zero-initialized table to store
   writer objects for the Fortran program.  */

#define VTK_XMLWRITERF_MAX 256
static vtkXMLWriterC* vtkXMLWriterF_Table[VTK_XMLWRITERF_MAX+1];

/**
 * \brief   vtkXMLWriterC_New Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_New(int* self)
{
  int i;

  /* Initialize result to failure.  */
  *self = 0;

  /* Search for a table entry to use for this object.  */
  for(i=1;i <= VTK_XMLWRITERF_MAX; ++i)
    {
      if(!vtkXMLWriterF_Table[i])
        {
          vtkXMLWriterF_Table[i] = vtkXMLWriterC_New();
          if(vtkXMLWriterF_Table[i])
            {
              *self = i;
            }
          return;
        }
    }
}

/**
 * \brief   vtkXMLWriterC_Delete Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_Delete(int* self)
{
  /* Check if the writer object exists.  */
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      /* Delete this writer object.  */
      vtkXMLWriterC_Delete(vtkXMLWriterF_Table[*self]);

      /* Erase the table entry.  */
      vtkXMLWriterF_Table[*self] = 0;
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_Delete called with invalid id %d.\n",
              *self);
    }

  /* The writer object no longer exists.  Destroy the id.  */
  *self = 0;
}

/**
 * \brief   vtkXMLWriterC_SetDataObjectType Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetDataObjectType(const int* self, const int* objType)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetDataObjectType(vtkXMLWriterF_Table[*self], *objType);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetDataObjectType called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetDataModeType Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetDataModeType(const int* self, const int* datamodetype)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetDataModeType(vtkXMLWriterF_Table[*self], *datamodetype);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetDataModeType called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetExtent Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetExtent(const int* self, int extent[6])
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetExtent(vtkXMLWriterF_Table[*self], extent);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetExtent called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetPoints Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetPoints(const int* self, const int* dataType, void* data,
   const vtkIdType* numPoints)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetPoints(vtkXMLWriterF_Table[*self], *dataType,
                              data, *numPoints);
    }
  else
    {
    fprintf(stderr,
            "vtkXMLWriterF_SetPoints_Float called with invalid id %d.\n",
            *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetOrigin Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetOrigin(const int* self, double origin[3])
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetOrigin(vtkXMLWriterF_Table[*self], origin);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetOrigin called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetSpacing Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetSpacing(const int* self, double spacing[3])
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetSpacing(vtkXMLWriterF_Table[*self], spacing);
    }
  else
    {
    fprintf(stderr,
            "vtkXMLWriterF_SetSpacing called with invalid id %d.\n",
            *self);
    }
}

/**
 * \brief   vtkXMLWriterC_SetCoordinates Fortran to C wrapper
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetCoordinates(const int* self, const int* axis,
   const int* dataType, void* data, const vtkIdType* numCoordinates)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetCoordinates(vtkXMLWriterF_Table[*self], *axis,
                                   *dataType, data, *numCoordinates);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetCoordinates called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetCellsWithType
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetCellsWithType(const int* self, const int* cellType,
   const vtkIdType* ncells, vtkIdType* cells, const vtkIdType* cellsSize)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetCellsWithType(vtkXMLWriterF_Table[*self], *cellType,
                                     *ncells, cells, *cellsSize);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetCellsWithType called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetCellsWithTypes
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetCellsWithTypes(const int* self, int* cellTypes,
   const vtkIdType* ncells, vtkIdType* cells, const vtkIdType* cellsSize)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetCellsWithTypes(vtkXMLWriterF_Table[*self], cellTypes,
                                      *ncells, cells, *cellsSize);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetCellsWithTypes called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetPointData
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetPointData(const int* self, const char* name,
  const int* dataType, void* data, const vtkIdType* numTuples,
  const int* numComponents, const char* role)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetPointData(vtkXMLWriterF_Table[*self], name, *dataType,
                                data, *numTuples, *numComponents, role);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetPointData called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetCellData
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetCellData(const int* self, const char* name,
  const int* dataType, void* data, const vtkIdType* numTuples,
  const int* numComponents, const char* role)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetCellData(vtkXMLWriterF_Table[*self], name, *dataType,
                                data, *numTuples, *numComponents, role);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetCellData called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetFileName
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetFileName(const int* self, const char* fileName)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetFileName(vtkXMLWriterF_Table[*self], fileName);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetFileName called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_Write
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_Write(const int* self, int* success)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      *success = vtkXMLWriterC_Write(vtkXMLWriterF_Table[*self]);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_Write called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_SetNumberOfTimeSteps
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_SetNumberOfTimeSteps(const int* self, const int* numTimeSteps)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_SetNumberOfTimeSteps(vtkXMLWriterF_Table[*self],
                                         *numTimeSteps);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_SetNumberOfTimeSteps called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_Start
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_Start(const int* self)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_Start(vtkXMLWriterF_Table[*self]);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_Start called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_WriteNextTimeStep
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_WriteNextTimeStep(const int* self, const double* timeValue)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_WriteNextTimeStep(vtkXMLWriterF_Table[*self], *timeValue);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_WriteNextTimeStep called with invalid id %d.\n",
              *self);
    }
}

/**
 * \brief   vtkXMLWriterF_Stop
 * \author  Ken Martin, Will Schroeder, Bill Lorensen
 */

void vtkXMLWriterF_Stop(const int* self)
{
  if(*self > 0 && *self <= VTK_XMLWRITERF_MAX && vtkXMLWriterF_Table[*self])
    {
      vtkXMLWriterC_Stop(vtkXMLWriterF_Table[*self]);
    }
  else
    {
      fprintf(stderr,
              "vtkXMLWriterF_Stop called with invalid id %d.\n",
              *self);
    }
}
