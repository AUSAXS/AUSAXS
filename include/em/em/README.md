This folder contains the electron-microscopy classes used to validate EM maps against SAXS data.

An EM map is read as a stack of 2D images (`Image.h`, `ImageStack.h`). By thresholding the map density, a dummy atomic structure is constructed and its scattering profile compared to the experimental curve. `ObjectBounds2D`/`ObjectBounds3D` track the occupied region of the map, and `manager/` holds the protein managers that build the dummy structure and drive the threshold search.
