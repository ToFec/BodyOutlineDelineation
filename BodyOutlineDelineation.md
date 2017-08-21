# Fully Automatic Body Outline Delineation on CT

This algorithm detects the body outline fully automatically for CT input images. 
At the moment the algorithm works well for abdomen and thorax CT images. For head and neck there are still some problems.
The algorithm reads and writes NRRD-images.

Call the algorithm as follows:

    ./BodyOutlineDelineation InputImageName.nrrd OutputImageName.nrrd
