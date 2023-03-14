# Toward mesh adaptation

[toc]

## Preliminaries

- Revise `refinement` markers
  - `100` regular active
  - `50` regular frozen
  - `1:6` face frozen
  - `7:18` edge frozen
  - `19:26` vertex frozen
- Parent level
  - `element%adapdation%refinement` is set to old refinement
  - `element%adapdation%sublevels` is set to new number of subdivisions  
  - `element%adaptation%mark` is set  to  new refinement
  - `ghost%adaptation%mark` is set to new refinement of master

## Child level partitioning

1. Check if ParMetis interface can be unified: yes :white_check_mark:
2. New type ChildDistributionMap_3D :white_check_mark:

## Child generation

1. Child TPs and IDs of parent face neighbors :white_check_mark:
2. Child TPs and IDs of parent edge neighbors :white_check_mark:
3. Child TPs and IDs of parent vertex neighbors :white_check_mark:
4. Child element data, face data, neighbors and geometry :white_check_mark:
5. ?

## Child transfer

## Assignment of retained data

## Child level completion

