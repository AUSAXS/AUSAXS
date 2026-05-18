This folder contains the file input/output layer.

Contents
- `File.h`, `Folder.h`, `ExistingFile.h` — path abstractions that distinguish files from folders and validate existence.
- `Reader.h`, `Writer.h` — the generic structure-file read/write interfaces.
- `pdb/` — readers and writers for the PDB structure format.

Additional structure and trajectory format support lives in the (private) `detail/` folder.
