## Stokes_GMS

This  example is proposed in Guermond et al. (2006) for exploring the effect of splitting error in projection schemes. It is generalized here to 3D by imposing periodicity in the z-direction.

To run the test with default parameters type

    mpirun -n 1 ../isp_flow__sdc_test

The solver parameters can be changed by editing the input file `isp_flow__sdc_test.prm`. 

To facilitate parallel runs you may increase the number of partitions per direction, `np(1:3)`, and adjust the number elements per direction in each partition, `ep(1:3)`. The number of required MPI processes equals the product of the `np` components. For example, setting `np = 2,2,1` results in four partitions and, hence, the test has to be started with 

    mpirun -n 4 ../isp_flow__sdc_test

### References

J. Guermond, P. Minev & J. Shen, _An overview of projection methods for incompressible flows_, Computer Methods Appl. Mechanics & Engineering, 195(44-47):6011–6045, 2006