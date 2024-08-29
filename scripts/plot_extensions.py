from plot_helper import * 

def load_file(file: str):
    def determine_type(line: str) -> PlotType:
        for t in PlotType:
            if line == t.value:
                return t
        return PlotType.Invalid

    data = []
    with open(file) as f:
        # keep reading until eof
        while line := f.readline():
            # determine the type of plot to make
            match determine_type(line.rstrip()):
                case PlotType.Dataset:
                    dataset: Dataset = read_dataset(f)
                    if dataset.data.size == 0:
                        if dataset.options.title != "":
                            title = dataset.options.title
                        elif dataset.options.legend != "":
                            title = dataset.options.legend
                        else:
                            title = "(unnammed)"
                        print(f"Skipping empty dataset \"{title}\"")
                        continue
                    data.append(dataset)

                case PlotType.Hline:
                    hline: Hline = read_hline(f)
                    data.append(hline)

                case PlotType.Vline:
                    vline: Vline = read_vline(f)
                    data.append(vline)

                case PlotType.Landscape:
                    dataset: Dataset = read_dataset(f)
                    data.append(dataset)

                case PlotType.Histogram:
                    dataset: Dataset = read_dataset(f)
                    data.append(dataset)

                case PlotType.Image:
                    x, y, z, _ = read_2dhist(f)
                    data.append((x, y, z))

                case PlotType.ImageAtoms:
                    print("ImageAtoms not implemented yet.")
                    exit(1)

                case _:
                    print("plot_file: Invalid plot type: " + line)
                    exit(1)

    return data