# Loci Radial Distribution Analysis
 
Fiji macro for Loci Radial Distribution Analysis

Developed for Richard Mann Lab, Columbia University. 2019
Similar to methods used in https://lens.elifesciences.org/28975/

*Requires 3D ImageJ Suite
 
For the detection of gene loci and subsequent analysis of protein distribution around this loci.
Loci are detected and radial distribution is analyzed in channels of interest and average distribution is generated (normalised to 1) and background subtracted (lowest value - outside of nucleus - assumed to be 0). Plot of intensity normalised 1 but also variance curves of distribution can be plotted as spheres in original image, allowing for spatial analysis.

Update May 16th 2019
Automatic detection uses thresholded region - then extracts local volume, performs maxima detection, then plots centroid back into the data
Added 2D manual selection - which uses MaxIP but searches through full 3D around loci XY location.
Modified 3D manual selection - uses local search rather than full image seach - much faster 1min vs 20-60min
Add random seeding for validation - this is measuring only Channel of Interest (COI), futher modification could be made to create seeds and measure COI and Loci channel, and create visualizations and snapshots area surrounding centroids. Future analysis could include clustering analysis


1. copy the .jar file into your ImageJ/Fiji plugins folder and restart ImageJ.
2. Run from Plugins > Cellular Imaging > Loci Radial Distribution Analysis
3. Select the folder containing raw data
4. Select the method for loci detection. an interactive mode is available for difficult datasets
5. Select the channel containing the loci, and the channel for analysis
6. Set the radius to be used for the radial distribution analysis (3um works well)

7. When using Interactive Loci Detection, the channel containing loci will open and you will need to select the loci, the multi point tool will be activated so you need to scroll through the z-stack and click on the loci. You don't have to be super precise I've written the pipeline to use these as starting points, it will search nearby for the exact 3D center of the Loci.

8. After running the analysis you will find a subfolder called "Analyzed". In this folder will be several files for each original input image. All the filenames will have the prefix of the input image.

If the SNR is high enough you should be able to use auto-detection in most cases. In cases where auto-detection is not appropriate a 3D manual selection and a 2D manual selection are available. The 2D manual selection allows you to select the loci from a Maximum Intensity Projection image, but will subsequently search for the loci in 3D space - this is a fast manual approach and has worked well in my hands on the test data. In my hands it takes ~3 minutes to process a 2GB image.

The output from the pipeline includes:
- 3 tables recording the raw radial distribution measurements (measured on the channel of interest (COI) at the loci site, the loci channel at the same location, and the COI at some randomly seeded sites within the nuclei). - Image_Name_Raw_Radial_Distribution_Ch_x.csv
- 3 tables recording the normalized radial distribution measurements (normalized to center intensity with the minimum intensity subtracted - from a value outside the nucleus) - Image_Name_Normalized Radial_Distribution_Ch_x.csv
- 3 custom summary tables including the values we discussed: center intensity, 250nm sphere mean intensity, and surrounding nucleus mean intensity, along with ratios of these. -	Image_Name_Summary_of_Loci_Measurements_Ch_x.csv
- A validation image which shows the loci channel raw data, and the centroids of the detected loci - useful to validate the performance of the loci detection - 	Image_Name_Loci_validation.tif
- A 3D image with spheres plotted at the loci locations with their intensity being the normalized COI intensity at center (e.g. 0.7 = 70) allows easy visualization of sites that correlated with the loci. -Image_Name_Loci_Intensity_COI_Ch_x.tif
- A 3D image with spheres in the loci locations with their intensity being the sphere to nucleus ratio (x100, so a ratio of 2.52 = 252) - Image_Name_Loci_250_Sphere_Nuclueus_Ratio.tif
- An image file containing all the detected loci, each slice is an image at the center of the loci and the surrounding area. It's a good image to use to visualize each loci and compare with the measurements. - Image_Name_Loci_Locations.tif
