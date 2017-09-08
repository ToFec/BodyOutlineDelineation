# Fully Automatic Body Outline Delineation on CT

This algorithm detects the body outline fully automatically for CT input images. 
At the moment the algorithm works well for abdomen and thorax CT images. For head and neck there can be some problems. Still needs to be evaluated.
The algorithm reads and writes NRRD-images. Optional it can write a meta information file with the computation time

Call the algorithm as follows:

    ./BodyOutlineDelineation InputImageName.nrrd OutputImageName.nrrd MetaInfFileName

The meta information file is optional. 
