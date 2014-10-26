rgc
===
These files implement the methods described in the paper
U. Sümbül et al., ``A genetic and computational approach to structurally classify neuronal types,'' Nature Communications, 2014.

These files are authored and maintained by Uygar Sümbül unless otherwise noted. 

rgcAnalyzer.m -- This file includes various Matlab functions that obtain a registered representation of a retinal ganglion cell starting
from its arbor trace and starburst amacrine surface annotations, performs hierarchical clustering, and the cutting of the hierarchical
clustering tree in the presence of partial cluster membership information. Authorship and redistribution information for the publicly
available functions used in this file is provided at the beginning of individual functions.

GPUrecon_bands.m, GPUreconv2.m, bandCNN.mat, currentnetworkv5_net249.mat -- These files implement the convolutional neural networks that
enhance the retinal ganglion cell arbors and starburst amacrine surfaces. They require CUDA-capable GPUs and associated drivers, and the
installation of cnpkg2, which is a publicly available package for creating 3-D convolutional networks and training them via the
backpropagation algorithm, by Srini Turaga: http://github.com/srinituraga/cnpkg

NNParams.mat, NNpass.m, processchat.m, processchatNN.m -- These files are authored by Sen Song, Tsinghua University. They detect the
starburst amacrine surfaces given an enhanced stack. They include a basic machine-learning based filtering architecture and a
maximum-surface detection routine.

---------------------- Instructions on rgcAnalyzer.m ----------------------

This file performs various tasks. The comments in the first function (rgcAnalyzer) explain these tasks. The file expects 5 inputs from
the user: arborFileName, OnSACFilename, OffSACFilename, voxelRes, conformalJump.

arborFileName: the path of the swc file that contains the trace of the neuron to be analyzed. The raw stack should be oriented so that
the outer retina is at the top. (i.e., The On starburst surface has smaller depth values than the Off starburst surface.) The files
used in the paper were generated using the Simple Neurite Tracer plug-in in FIJI.

OnSACFilename: the path of the txt file that contains the annotations for the On starburst surface. (See below.)

OffSACFilename: the path of the txt file that contains the annotations for the Off starburst surface. (See below.)

voxelRes: physical xyz resolution of a voxel in the image stack in micrometers (1x3 array). The value used in the paper is [0.4 0.4 0.5].

conformalJump: (positive integer, typically 1, 2, 3) the length of the edge of the grid used in conformal mapping in terms of voxels.
"1" performs exact calculation, however it can be slow. The interpolation introduced by using "2" or "3" is usually minimal, and the
mapping runs faster.

---------------------- Instructions for generating surface annotation files in FIJI ----------------------

- Open the three-dimensional image stack containing the starburst cell staining in FIJI.
- Image>Stacks>Reslice. From the menu, choose "1 pixel output spacing", "start at top", "avoid interpolation."
- Choose the point tool from the main toolbar. Double click on the point tool to configure: auto-measure, label points.
- Start with the annotation of the On starburst layer - annotate only the top layer. Pick a slice and click on 5-10 points on the On
layer in that slice. This should open a Results window indicating the X,Y, and Slice positions. (Since the stack is resliced, Y is
actually Z. This is automatically handled in the code.)
- Make sure that the point annotation tool records pixel values and not physical values.
- Advance a few tens of slices and repeat the previous step. Continue until you reach the end of the stack.
- Save the On starburst annotations as a txt file from File> Save as under the Results window menu.
- Repeat the previous four steps for the Off starburst layer.
