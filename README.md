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
like to contribute to this project, please feel free to open an issue or submit a pull request.

For questions or inquiries, please contact:
- Dr. Shuo Pang: pang@ucf.edu
- Xichen Shan: xcshan@ucf.edu