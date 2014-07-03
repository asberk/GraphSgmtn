#README
-----

There are three main files in this project. 
* `PredImgSgmtn.m`
* `PredUSHouseVoting.m`
* `NonlocalMeans.m`

This project also comes with a few helper functions for debugging,
ease of data manipulation.
* `getImage.m`
* `imgnmz.m`
* `reOrder.m`

Lastly, this project has several supporting functions on which the
main code depends. 
* `gaussNorm.m`
* `getFeatureArr.m`
* `getGraphLaplacian.m`
* `getImage.m`
* `getNormalizedSampleWeights.m`
* `getWeights.m`
* `getWeightsXY.m`
* `imgnmz.m`
* `widen.m`

`PredImgSgmtn.m` reads in an image and uses the Nystr\"om
extension to approximately solve for the adjacency matrix for the
given input image. The adjacency matrix is a weight matrix,
representing the edge weights between each pixel, for which an
undirected graph has been generated, and in which each vertex has
associated to it a feature vector/matrix, which is given by the
NxN-neighbourhood of pixels about that pixel [vertex]. The
eigenvectors of the adjacency matrix (approximated by the Nystrom
extension) are used to quickly solve a nonlinear diffuse interface
problem, posed using a Ginzburg-Landau-type functional with a
fidelity term. The solution of this Ginzburg-Landau-type
functional represents a two-phase solution, to which thresholding
techniques may be applied to yield a final segmentation. 

`PredUSHouseVoting.m` follows a similar structure; however, the
dataset is small enough to compute the full weight matrix, and
obtain its true eigenvalues, without additional approximation
methods. Given a small set of samples (5), this dataset predicts
the party affiliation of the remaining (430) samples with
impeccable accuracy --- again, by solving a nonlinear PDE designed
to represent graph cutting. 

All files come with help documentation. Please see the Matlab help
documentation on any file `myFile.m` by typing `help myFile` in
the Matlab command prompt. 



