This folder contains the dataset classes that hold tabular data — most importantly experimental scattering curves (q, I, and uncertainties).

Contents
- `Dataset.h`, `Dataset2D.h`, `SimpleDataset.h` — the core dataset types, ranging from a bare 2D table to a scattering curve with errors.
- `PointSet.h`, `Multiset.h` — collections of points and of datasets.
- `NamedDataset.h`, `NamedWrapper.h` — attach labels to datasets and their columns.
- `DatasetFactory.h` — constructs datasets by reading data files from disk.
