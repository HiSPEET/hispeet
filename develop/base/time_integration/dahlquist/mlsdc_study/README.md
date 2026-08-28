# MLSDC studies for the Dahlquist problem

## Interactive shell

- Select example

  - Euler-based method

    ```bash
    cd mlsdc_eu_eu
    ```

  - 1-stage implicit streamline diffusion

    ```bash
    cd mlsdc_s1_s1
    ```

  - 2-stage implicit streamline diffusion

    ```bash
    cd mlsdc_s2_s2
    ```

- Conduct study

  - Run in shell

    ```bash
    bash -l mlsdc_study.sh
    ```

  - Run in batch mode using `slurm` 

    ```bash
    # intialize environment, e.g, on barnard.hpc.tu-dresden.de
    ml barnard/gnu/12
    sbatch mlsdc.slurm
    ```

  - Be patient: the studies run for several minutes
