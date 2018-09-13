/**
 * \file       exodus_wrapper.c
 * \brief      Wrappers for accessing selected ExodusII functions from Fortran
 * \author     Joerg Stiller
 * \date       2016/05/21
 * \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
 *
 * \details
 * ...
 ******************************************************************************/

/* Exodus C interface */
#include "exodusII.h"

/**
 * ex_open
 */

int exf_open (const char* path, int mode, int* comp_ws, int* io_ws,
              float* version)
{
  return ex_open (path, mode, comp_ws, io_ws, version);
}

/**
 * ex_close
 */

void exf_close (int exoid, int* err)
{
  *err = ex_close (exoid);
}

/**
 * ex_get_init
 */

void exf_get_init (int exoid, char* title, int* num_dim, int* num_nodes,
                   int* num_elem, int* num_elem_blk, int* num_node_sets,
                   int* num_side_sets, int* err)
{
  *err = ex_get_init (exoid, title, num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);
}

/**
 * ex_get_coord -- generic version
 */

void exf_get_coord (int exoid, void* x_coor, void* y_coor, void* z_coor,
                    int* err)
{
  *err = ex_get_coord (exoid, x_coor, y_coor, z_coor);
}

/**
 * ex_get_elem_block
 */

void exf_get_elem_block (int exoid, int elem_blk_id, char* elem_type,
                         int* num_elem_this_blk, int* num_nodes_per_elem,
                         int* num_attr, int* err)
{
  *err = ex_get_elem_block (exoid, elem_blk_id, elem_type, num_elem_this_blk,
                            num_nodes_per_elem, num_attr);
}

/**
 * ex_get_elem_blk_ids
 */

void exf_get_elem_blk_ids (int exoid, int* elem_blk_ids, int* err)
{
  *err = ex_get_elem_blk_ids (exoid, elem_blk_ids);
}

/**
 * ex_get_elem_conn
 */

void exf_get_elem_conn (int exoid, int elem_blk_id, int* connect, int* err)
{
  *err = ex_get_elem_conn (exoid, elem_blk_id, connect);
}

/**
 * ex_get_name
 */

void exf_get_name (int exoid, int obj_type, int id, char* name, int* err)
{
  *err = ex_get_name (exoid, obj_type, id, name);
}

/**
 * ex_get_side_set_ids
 */

void exf_get_side_set_ids (int exoid, int* side_set_ids, int* err)
{
  *err = ex_get_side_set_ids (exoid, side_set_ids);
}

/**
 * ex_get_side_set_param
 */

void exf_get_side_set_param (int exoid, int side_set_id, int* num_side_in_set,
                             int* num_dist_fact_in_set, int* err)
{
  *err = ex_get_side_set_param (exoid, side_set_id, num_side_in_set,
                                num_dist_fact_in_set);
}

/**
 * exf_get_side_set
 */

void exf_get_side_set (int exoid, int side_set_id, int* side_set_elem_list,
                             int* side_set_side_list, int* err)
{
  *err = ex_get_side_set (exoid, side_set_id, side_set_elem_list,
                          side_set_side_list);
}


