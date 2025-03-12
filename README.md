# MATLAB FDFD Photonic Simulation Toolbox
Photonic Finite-Difference Frequency-Domain (FDFD) Toolbox is a MATLAB Community Toolbox project, aiming 
to provide intuitive framework for photonic devices simulation using MATLAB.

> **Note:** The toolbox is currently in its early development stage. New features for high-dimensional 
simulation, advanced analysis, and optimization will be coming in next release!

## Key Features
- **Customizable Geometry**: Simple geometry and mesh creation.
- **Object-Oriented Implementation**: Leverage modular MATLAB Object-Oriented Programming (OOP) capability
for expandability and scalability.

## Prerequires
- **MATLAB**: Version R2014b or later

## Toolbox Installation
1. Clone/download this repository.
2. Start MATLAB.
3. Install the Toolbox by double-clicking on the `FDFDPhotonicToolbox.mltbx` in MATLAB.

## Examples 
Under the `examples/` directory, live scripts demonstrate the working flow of the Toolbox:

| Dimension | Structure | Description | File Name | 
|---|---|---|---|
| 2D | Directional coupler | 2D Scattering | `directional_coupler.mlx` |
| 2D | Grating coupler |2D Scattering | `grating_coupler.mlx` |
| 2D | 1x2 Inverse-designed mode division multiplexer |2D Scattering | `mode_mux_2x1.mlx` |
| 2D | Transmission line | 2D Eigenmode (metal) | `transmission_line.mlx` |

\* For the example file `mode_mux_2x1.mlx`, GDSII Toolbox is required (see section [GDSII Toolbox](#gdsii-toolbox)).



## Quick Start  
For quick start, please checkout the example models under `.\examples` directory. Here we list some frequently used methods:

| Top Class  | Method Syntax                                   | Description |
|------------|----------------------------------------------|-------------|
| **Layout** | `setGeometry('Rectangle', L, W, 'option')` | Create a rectangel geometry with length `L` and width `W`. |
|            | `setGeometry('Disk', R, C, 'option')`      | Create a disk geometry with radius `R` at center coordinates `C`. |
|            | `setGeometry('Polygon', X, Y)`            | Create a plygon geometry with at a list of vertices at coordinates array `X` and `Y`. |
|            | `addPort(dir, position, range)`          | Add a port with its norm in `dir` direction, located at `position`, with spatial extent defined in `range`. |
| **Material** | `setMaterial(eps_r)`                      | Assign the material with relative permittivity of `eps_r`. |
| **Grid**     | `setMesh([Nx, Ny], [dx, dy])`            | Configure the simulation mesh grid with dimensions `Nx*Ny` and pixel size `dx*dy`. |
| **Source**   | `addSource(type, 'option')`             | Add an excitation source to the simulation  |



## GDSII Toolbox
GDSII is the industry standard for layouts, we have added the support for `.gds` files in our toolbox for user's convenience. 
We utilize the [Octave / MATLAB Toolbox for GDSII Stream Format](https://github.com/ulfgri/gdsii-toolbox) by Ulf Griesmann.

 To use the GDSII toolbox, please compile the mex functions in MATLAB on Windows by changing to the `./gdsii-toolbox` directory and running

```
>> makemex
```
To add the GDSII Toolbox to MATLAB path
```
>> addpath(genpath(<path/to/gdsii-toolbox>));
>> savepath;
```

## Feedback & Contributions
We appreciate your feedback and contributions! If you encounter issues, have suggestions, and would 
like to contribute to this project, please feel free to open an issue or contact us.

For questions or inquiries, please contact:
- Dr. Shuo Pang: pang@ucf.edu
- Xichen Shan: xcshan@ucf.edu


## References
[1] Rumpf, R. C. (2022). 
Electromagnetic and Photonic Simulation for the Beginner: Finite-Difference Frequency-Domain in MATLAB®. Artech House.

[2] Rumpf, R. C. (2012). 
Simple implementation of arbitrarily shaped total-field/scattered-field regions in finite-difference frequency-domain. Progress In Electromagnetics Research B, 36, 221-248.

[3] Christ, A., & Hartnagel, H. L. (1987). Three-dimensional finite-difference method for the analysis of microwave-device embedding.
 IEEE Transactions on Microwave Theory and Techniques, 35(8), 688-696.

[4] Taflove, A., Hagness, S. C., & Piket-May, M. (2005). Computational electromagnetics: the finite-difference time-domain method. The Electrical Engineering Handbook, 3(629-670), 15.

[5] Berenger, J. P. (1994). A perfectly matched layer for the absorption of electromagnetic waves. Journal of computational physics, 114(2), 185-200.

[6] W. Shin and S. Fan, “Choice of the perfectly matched layer boundary condition for frequency-domain Maxwell's equations solvers,” 
Journal of Computational Physics, vol. 231, pp. 3406–3431 (2012)

[7] W. Shin, MaxwellFDFD Webpage, 2015. https://github.com/wsshin/maxwellfdfd.