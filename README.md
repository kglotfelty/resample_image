# `resample_image`

Synopsis: _reproject by resampling_

## Concept: reprojection

Reprojection is described here as the transformation of an image 
with one Coordinate System (CS) to a different CS.  Often in astronomical
images we talk about the World Coordinate System (WCS), that is the
transformation from pixel(x,y) to celestial coordinates: right assenstion and
delication.

    $$$
    (X,Y) \rightarrow^{WCS} \alpha,\delta
    $$$

