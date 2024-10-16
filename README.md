# Slabify
A CLI tool to automatically segment the lamella slab in cryo-ET volumes

[![DOI](https://zenodo.org/badge/848790613.svg)](https://doi.org/10.5281/zenodo.13941082)
![slabify](https://github.com/user-attachments/assets/7c30cf40-76be-4293-ab65-dd5a651ced6b)


# Overview
Slabify is a Python command-line tool to automatically segment the lamella slab from cellular cryo-electron tomography (cryo-ET) volumes. The program analyzes the local variance around random points inside the tomogram to find where the "interesting" density is. There are three modes of operation:
1. Find the lamella boundaries by fitting two planes to the top and bottom sides of the slab iteratively (default).
2. Fit a single plane through the center of the lamella, then expand a slab mask of the given `--thickness` in Z. This mode is enabled by the `--simple` flag, and tends to work better in "difficult" cases.
3. If all else fails, you can still manually define the lamella boundaries by clicking a few points (12, to be precise) in IMOD and have `slabify` fit two planes defining the top and bottom sides of the slab. See [instructions](https://github.com/CellArchLab/slabify-et/wiki/How-to-manually-create-a-lamella-boundary-mask-using-IMOD) for details.

The first two modes employ robust fitting through the [RANSAC](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RANSACRegressor.html) algorithm to prevent outliers from perturbing the boundary estimation too much.

## Applications
A lamella slab mask can be useful in different scenarios when processing cellular cryo-ET data, for example:

* Removing obvious false positives outside the lamella in particle picking (e.g. by template matching)
* Speeding up template matching by avoiding searches in the void, as implemented in [STOPGAP](https://github.com/wan-lab-vanderbilt/STOPGAP), [GAPSTOP](https://gitlab.mpcdf.mpg.de/bturo/gapstop_tm) and [pytom-match-pick](https://github.com/SBC-Utrecht/pytom-match-pick)
* Clipping segmentations
* Selecting only "interesting" subtomograms for training denoising networks, as in [IsoNet](https://github.com/IsoNet-cryoET/IsoNet) and [DeepDeWedge](https://github.com/MLI-lab/DeepDeWedge)
* Assessing particle distribution within the lamella

# Installation

1. Create a conda environment (or use an existing one if you prefer):
```bash
conda create -n slabify python=3.12
```
2. Activate the environment:
```bash
conda activate slabify
```
3. Clone the repo and install Slabify using `pip`:
```bash
git clone https://github.com/CellArchLab/slabify-et.git
cd slabify-et
pip install .
```
Now, you should be able to run the `slabify` command.

## Optional dependencies
* [IMOD](https://bio3d.colorado.edu/imod/): if you want to be able to create slab masks from IMOD models in `.mod` format directly, you should have IMOD installed and activated in your session.

# Usage
Using Slabify is straightforward:
```bash
slabify --input tomogram.mrc --output tomogram_slab_mask.mrc
```
You may want to turn a few knobs to fine tune your results of course. Type `slabify --help` to see all options. See also the [FAQ](https://github.com/CellArchLab/slabify#faq) below.

# FAQ
### What tomogram dimensions are good for Slabify?
The following tomogram sizes have been tested:

* 1024 x 1024 x 512, 7.84 Å/px (bin4, Falcon4i detector)
* 1024 x 1440 x 256, 10.74 Å/px (bin4, K3 detector)

Anything in this range of tomogram dimensions and pixel sizes should work well. The `--boxsize` and `--n-samples` options might need to be adjusted if you have something very different.

### What kind of filtering should I apply to my tomograms?
While Slabify should work on any tomogram, high-contrast and "clean" tomos tend to give the best results. These can be obtained by applying a denoising tool such as [cryoCARE](https://github.com/juglab/cryoCARE_pip), [IsoNet](https://github.com/IsoNet-cryoET/IsoNet) or [DeepDeWedge](https://github.com/MLI-lab/DeepDeWedge) to your data before running Slabify. For a quicker alternative, the deconvolution filter from [Warp]((https://doi.org/10.1038/s41592-019-0580-y)) also works well for our purposes. See [here](https://github.com/CellArchLab/slabify-et/wiki/How-to-deconvolve-a-tomogram-using-IMOD) how to apply it using IMOD.

### The estimated slab is too thin, how can I improve it?
It is known that Slabify can be quite conservative in its estimation of the lamella boundaries. There are a few tips you can try to get thicker slab masks:
* The easiest way is to use the `--offset` option to arbitrarily grow your slab mask in the Z direction. Note that a negative value can be provided to make the slab *thinner*!
* You can increase the number of `--iterations`. 3 to 5 iterations are usually good.
* You can slightly decrease the `--percentile` of highest variance points in order to capture more information, say from 95 to 94.
* Tweaking the `--n-samples` and `--boxsize` values can be beneficial in some cases (needs more testing).
* Finally, the `--simple` mode with an arbitrary slab `--thickness` is generally a safe option.

# Acknowledgments
The manual boundary masking is based on code from Will Wan for [STOPGAP](https://github.com/wan-lab-vanderbilt/STOPGAP). Thanks to Caitie McCafferty and Philippe Van der Stappen for discussions and feature suggestions.
