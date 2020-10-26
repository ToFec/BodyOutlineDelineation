# Fully Automatic Body Outline Delineation on CT

This algorithm detects the body outline fully automatically for CT input images. 
At the moment the algorithm works well for abdomen and thorax CT images. For head and neck there can be some problems. Still needs to be evaluated.
The algorithm reads and writes NRRD-images. Optional it can write a meta information file with the computation time

Call the algorithm as follows:

    ./BodyOutlineDelineation InputImageName.nrrd OutputImageName.nrrd MetaInfFileName

The meta information file is optional.

In case you find this programm useful and use it for your work please cite: Fechter, Tobias & Dolz, Jose & Nestle, Ursula & Baltas, Dimos. (2017). EP-1417: Clinical evaluation of a fully automatic body delineation algorithm for radiotherapy. Radiotherapy and Oncology. 123. S757-S758. 10.1016/S0167-8140(17)31852-2. 
