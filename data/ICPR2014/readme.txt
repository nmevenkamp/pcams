
/------------------------------------------------------------------------------/
/ The Prague Texture Segmentation Data Generator and Benchmark readme file ;-) /
/------------------------------------------------------------------------------/

 Downloaded zip file includes this readme.txt, data.xml, edge maps, ground truth
and texture mosaics files in the PNG format. Generated data consist of several
sets. Each set is defined by its parameters determining the number of masks,
the number of texture mosaics per masks and the mask and textures properties.

 The names of image files are composed of a prefix which determines the type
of the file ('gt'=ground truth, 'mask'=edge file, 'tm'=texture mosaic,
'seg'=segmentation). The number of the set, the number of mask, the number
of mosaic (in case of 'tm' or 'seg') and the PNG extension are appended
to this prefix. For example the first texture mosaic name is 'tm1_1_1.png',
and the corresponding ground truth, edge and segmented images are 'gt1_1.png',
'mask1_1.png' and 'seg1_1_1.png', etc.

 Training texture set is included for classification (supervised segmentation).
For each texture mosaic there is a list of training textures (tags 'Train#.*')
associated with texture classes in the given mosaic. Class indeces (#) should
be used as values in the resulting segmentation image.

 You can upload your results for the online results evaluation. It is necessary
to upload zipped 'data.xml' and 'seg*.png' files (one 'seg' file for each 'tm'
file). The web interface url is http://mosaic.utia.cas.cz. You can evaluate
your results as a guest but you have to be a registered user if you want
to store and compare your results with others.

 If you have questions or comments please write an email to xaos@utia.cas.cz.


/-------- April 2010 --------/
Multispectral satellite texture set has been added. Generated texture mosaics
are in the PNG format with the postfix denoting the spectral band number.
More details can be found in the http://mosaic.utia.cas.cz/other/ali.txt file.



/-- data.xml description ------------------------------------------------------/

 The 'data.xml' file stores dataset parameters in a xml-like format. It has
'section' and 'var' tags. The structure of this file is described as follows:

'main':
  - 'Label'        (name of the benchmark)
  - 'NumOfSets'    (number of sets)
  + 'Set_*'        (set parameters)
  + 'params'       (parameters used for evaluation).

'Set_*':
  - 'Masks'        (number of masks)
  - 'TMsPerMask'   (number of texture mosaics per mask)
  - 'Note'         (set remark)
  + 'Mask_*'       (mask parameters)
  + 'GTMap_*'      (ground truth parameters)
  + 'TMMap_*_*'    (texture mosaic parameters)

'Mask_*':
  - 'File'         (mask filename)
  - 'Segments'     (number of segments)
  - 'Seed'         (random seed)
  - 'PointsX[]'    (x coordinates of the points generating mask)
  - 'PointsY[]'    (y coordinates of the points generating mask)
  - 'Rows'         (number of image rows)
  - 'Columns'      (number of image columns)
  - 'Type'         (mask type 0=straight 1=broken lines 2=splines)

'GTMap_*':
  - 'File'         (ground truth filename)
  - 'MaskRef'      (reference to Mask_*)

'TMMap_*_*':
  - 'File'         (texture mosaic filename)
  - 'BW'           (monospectral textures ? [0/1])
  - 'Key'          (invariant textures flags)
  - 'NN'           (number of invariant subregions)
  - 'GTMapRef'     (reference to GTMap_*)
  - 'SegmFile'     (filename of the segmentation result image)
  - 'Seed'         (random seed)
  - 'Textures[]'   (texture indices)
  - 'Train#.*'     (train texture(s) for class # - only for classification)


/-- data.xml example ----------------------------------------------------------/

<main>
  <var name="Label">"Benchmark dataset - Colour [normal] CL"</var>
  <var name="NumOfSets">20</var>
  <section name="Set_1">
    <var name="Masks">1</var>
    <var name="TMsPerMask">1</var>
    <var name="Note">"#3 - all"</var>
    <section name="Mask_1">
      <var name="File">"mask1_1.png"</var>
      <var name="Segments">3</var>
      <var name="Seed">1</var>
      <var name="PointsX[3]">216.000000 359.000000 107.000000</var>
      <var name="PointsY[3]">396.000000 236.000000 229.000000</var>
      <var name="Rows">512</var>
      <var name="Columns">512</var>
      <var name="Type">1</var>
    </section>
    <section name="GTMap_1">
      <var name="File">"gt1_1.png"</var>
      <var name="MaskRef">"Mask_1"</var>
    </section>
    <section name="TMMap_1_1">
      <var name="File">"tm1_1_1.png"</var>
      <var name="BW">0</var>
      <var name="Key">0</var>
      <var name="NN">1</var>
      <var name="GTMapRef">"GTMap_1"</var>
      <var name="SegmFile">"seg1_1_1.png"</var>
      <var name="Seed">1</var>
      <var name="Textures[3]">"99.A" "46.A" "52.A"</var>
      <var name="Train0.A" rows="512" cols="512">"textile/cloth29_t0.png"</var>
      <var name="Train1.A" rows="512" cols="512">"man-made/mat2_t0.png"</var>
      <var name="Train2.A" rows="512" cols="512">"man-made/roofTiles2_t0.png"</var>
    </section>
  </section>
  <section name="Set_2">
    .
    . And so on ...
    .
  </section>
  <section name="params">
    .
    . Parameters for evaluation purposes. Don't change them.
    .
  </section>
</main>

/------------------------------------------------------------------------------/
