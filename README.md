# VOICE (Vlasov Open boundary Ion Coupling to Electrons)

This repository contains 2D plasma simulations (1 spatial and 1 velocity dimension). The simulations are built using the gyselalibxx library [^gyselalibxx]. The simulations implemented are:
- Sheath simulation (see [^Bourne2023] for more details)
- Neutrals simulation

## Installing

## Building

```
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=gyselalibxx/toolchains/<MACHINE>/toolchain.cmake
cmake --build build -j 4
cmake --install build
```

## Running

```
```

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

***

[^gyselalibxx]: https://github.com/gyselax/gyselalibxx

[^Bourne2023]: Bourne, E., Munschy, Y., Grandgirard, V., Mehrenberger, M., & Ghendrih, P. (2023). Non-uniform splines for semi-lagrangian kinetic simulations of the plasma sheath. Journal of Computational Physics, 488, 112229. https://doi.org/10.1016/j.jcp.2023.112229 
