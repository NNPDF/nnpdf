# How to run a fit

By running a fit one performs a PDF determination using the NNPDF methodology. In a schematic way, one will generate a PDF from a grid on x with a neural network, then convolute the result with an FK table containing the partonic interaction, and finally compare with data to optimize a given cost function. For that one has to install the nnpdf code, then prepare a runcard and run the fit itself.

0. Install the code (follow NNPDF guide in this same repository)

1. Create a runcard by taking as template one of the files in `<profile_prefix>/config`. The `<profile_prefix>` path is by default `<install prefix>/share/NNPDF` for source installation, while `<conda root>/share/NNPDF` for conda installation.

2. Prepare the fit: use the command`vp-setupfit <runcard>.yml` to generate a `<runcard_folder>` folder in the current directory with a copy of the original YAML runcard. The required resources (such as the theory ID and the PDF) will be downloaded automatically. Alternatively they can be obtained
with the `vp-get` tool.

3. The `nnfit` program takes a `<runcard_folder>` as input, e.g.  ```nnfit <replica_number> <runcard_folder> ``` where replica_number goes from 1 to n. (You can change the o)

4. Wait until you have fit results, then use `postfit <number_of_replicas> <runcard_folder>` to finalize the PDF set by applying post selection criteria. This will produce a set of `<number_of_replicas>+1` replicas.

5. Upload the results using `vp-upload --fit <runcard_folder>`. Then, to analyze results with `validphys`, see the [vp-guide](https://data.nnpdf.science/validphys-docs/guide.html#development-installs). Consider running the `vp-comparefits -i` command.