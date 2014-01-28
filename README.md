rgc
===
These files implement the methods described in the paper
U. Sümbül et al., ``A genetic and computational approach to structural classification of neuronal types,'' Nature Communications, 2014.

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
