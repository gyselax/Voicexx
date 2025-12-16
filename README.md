# VOICE (Vlasov Open boundary Ion Coupling to Electrons)

This repository contains 2D plasma simulations (1 spatial and 1 velocity dimension). The simulations are built using the gyselalibxx library [^gyselalibxx]. The simulations implemented are:
- Sheath simulation (see [^Bourne2023] for more details)
- Neutrals simulation

## Installing

First, clone this repository and run the following command to initialise the submodules.
```
git submodule update --recursive --init
```
Then we need to source the correct environment.
For a list of available environments, see the [toolchains](./gyselalibxx/toolchains) folder,
and particularly the following [documentation](https://gyselax.github.io/gyselalibxx/toolchains/index.html#available_systems).

In the particular case of persee, run
```
. gyselalibxx/toolchains/persee/v100/environment.sh
```

## Building

For building, create a build folder and run
```
cmake -S <VOICE> -B build -DCMAKE_TOOLCHAIN_FILE=gyselalibxx/toolchains/<MACHINE>/toolchain.cmake
cd build && make -j 4
```

## Running

After compilation, the executables will be available in the [bin](./bin) folder.
To execute them, enter the following commands.
```
./sheath_xnonperiod_vx --dump-config param.yaml
./sheath_xnonperiod_vx param.yaml
```

Post-process functions are available in [post-process](./post-process/).
For documentation see this [readme](./post-process/README.md).

## Non-uniform grids

A tool for generating non-uniform grids is available in [non-unif-grids](./pre-process/suggested_points_refinement.py).
An example of usage is below.
```
python3 <VOICE>/pre-process/suggested_points_refinement.py grid.h5 \
    --name breakpoints_x --xmin 0 --xmax <grid_length> \
    --edge-domains 0 $p1 $p2 $p3 $p4 <grid_length> \
    --ncells 128 1024 2048 1024 128
python3 <VOICE>/pre-process/suggested_points_refinement.py maill_vx.h5 \
    --name breakpoints_vx --xmin -6 --xmax 6 --edge-domains -6 6 --ncells 256
h5copy -i maill_vx.h5 -o grid.h5 -s breakpoints_vx -d breakpoints_vx
rm maill_vx.h5
```

***

[^gyselalibxx]: https://github.com/gyselax/gyselalibxx

[^Bourne2023]: Bourne, E., Munschy, Y., Grandgirard, V., Mehrenberger, M., & Ghendrih, P. (2023). Non-uniform splines for semi-lagrangian kinetic simulations of the plasma sheath. Journal of Computational Physics, 488, 112229. https://doi.org/10.1016/j.jcp.2023.112229 
