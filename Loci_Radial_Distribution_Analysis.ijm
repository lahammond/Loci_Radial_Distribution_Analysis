// Cellular Imaging : Loci Radial Distribution Analysis
//
// Author: 	Luke Hammond
// Cellular Imaging | Zuckerman Institute, Columbia University
// Version: 0.9.12
// Date:	1st December 2018
//
// Developed for Rebecca Delker / Richard Mann lab Zuckerman Institute, Columbia University
//
// For the detection of gene loci and subsequent analysis of protein distribution around this loci.
// loci are detected and radial distribution is analyzed in channels of interest
// and average distribution is generated (normalised to 1) and background subtracted (lowest value - outside of nucleus - assumed to be 0)
// 		needs to plot intensity normalised 1 but also variance
// curves of distribution can be plotted as spheres in original image, allowing for spatial analysis

// Update May 16th 2019
// Automatic detection uses thresholded region - then extracts local volume, performs maxim detection, then plots centroid back into the data
// Added 2D manual selection - which uses MaxIP but searches through full 3D around loci XY location.
// Modified 3D manual selection - uses local search rather than full image seach - much faster 1min vs 20-60min
// Add random seeding for validation - this is measuring only Channel of Interest (COI), futher modification could be
// made to create seeds and measure COI and Loci channel, and create visualizations and snapshots area surrounding centroids


// Future analysis could include clustering analysis of loci with +ve correlation with ubx



// Initialization
requires("1.51w");



#@ File[] listOfPaths(label="select files or folders", style="both")
#@ String(label="Method for loci detection:", choices={"Automatic Detection", "2D Manual Loci Selection", "3D Manual Loci Selection"}, value = "Manual Loci Selection", style="listBox", description="") LociDetType

#@ String(label="Select loci channel for detection:", choices={"1", "2", "3", "4", "0"}, style="radioButtonHorizontal", value = "1", description="") LociCh

#@ String(label="Select channel for radial distribution analysis:", choices={"1", "2", "3", "4", "0"}, style="radioButtonHorizontal", value = "0", description="") DistCh1
#@ BigDecimal(label="Select radius for radial distribution (um):", value = 3, style="spinner") RDRad
#@ BigDecimal(label="Select minimum intensity for randomly seeded faux loci (for validation):", value = 800, style="spinner") SeedMinInt



run("Options...", "iterations=3 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=yellow");
run("Point Tool...", "type=Hybrid color=Yellow size=Medium label show counter=0");
run("Clear Results");

setBatchMode(true);

Isotropic = 0;
//CloseOpenWindows();


// Preparation

for (FolderNum=0; FolderNum<listOfPaths.length; FolderNum++) {
	input=listOfPaths[FolderNum];
	if (File.exists(input)) {
        if (File.isDirectory(input)) {

			// Create folders
			File.mkdir(input + "/Analyzed");
			ChOut = input + "/Analyzed/";
			LociOut = input + "/Temp/Loci/";
			MeasuredOut = input + "/Temp/Measured/";


			// get files
			files = getFileList(input);	
			files = ImageFilesOnlyArray(files);		
			
			run("Collect Garbage");

			//iterate over all files
			
			for(i=0; i<files.length; i++) {				
				starttime = getTime();
				closewindow("Log");
				File.mkdir(input + "/Temp");
				TempOut = input + "/Temp/";
				run("Clear Results");
				print("\\Clear");
		
				print(" Processing folder "+(FolderNum+1)+" of "+listOfPaths.length+" folders selected for processing." );
	       		print(" Processing folder "+(FolderNum+1)+": " + input + " ");
				
				
				
				image = files[i];	
				print("\\Update3: Processing image " + (i+1) +" of " + files.length +". Image:"+image+"");
				run("Bio-Formats Importer", "open=[" + input + "/" + image + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
				getVoxelSize(W, H, depth, unit);
				getDimensions(width, height, ChNum, slices, frames);

				File.mkdir(input + "/Temp/Loci");
				File.mkdir(input + "/Temp/Measured");
				
				// set Radial distribution step sizes
				RDRadSteps = parseInt(((RDRad*1000)/(W*1000)));
				// Calculate steps to average to measure mean int for 250nm volume
				sphere250 = parseInt(((250)/(W*1000)));
				// calculate rectangle area to measure approximately 
				box4um = parseInt(4/W);
				// Create box size for extraction - Extraction box MUST be even or relocalizing position will be inaccurate
				CentroidSearchSize = parseInt(1.3/W);
				if ( CentroidSearchSize % 2 > 0 ) {
					CentroidSearchSize ++;
				}
				
				if (Isotropic == 1) {
					CentroidSearchSizeZ = parseInt(1.3/W);
				} else {
					CentroidSearchSizeZ = parseInt(1.3/depth);
				}

				if ( CentroidSearchSizeZ % 2 > 0 ) {
					CentroidSearchSizeZ ++;
				}
				
				ExtractROIsize = parseInt(3.5/W);
								
				// calculate border region (2um on XY and ~1.5um Z)
				XY_edge = parseInt(2.1/W);
				//Z_edge = parseInt(1.5/depth); - since we are making isotropic in later step Z_edge can = XY_edge 
				if (Isotropic == 1) {
					Z_edge = XY_edge;
				} else {
					Z_edge = parseInt(2.1/depth);
				}
				//Extraction Region for 3D radial failures
				RadExtract = parseInt(2/W);
				reslicefactor = depth/W;
								
				
				print("Radius for radial distrubtion is set to "+RDRad+" micron.");
				print("Lateral image resolution is "+W+" microns. Axial resolution is "+depth+" microns.");
				print(RDRadSteps+"steps required for radial distribution radius.");
				print("Ignoring foci within "+XY_edge+" isotropic voxels of edge of 3d volume (~2.1 micron).");
			
				imagename = short_title(image);
				
				rename("Raw");
				run("Split Channels");

				selectWindow("C" + DistCh1 + "-Raw");
				rename("RDCh1");
				RDCh1ID = getImageID();
				selectWindow("C" + LociCh + "-Raw");
				rename("Loci");
				// Clean up other channels
				for (cs = 0; cs < 5; cs++) {
					closewindow("C" + cs + "-Raw");
				}
				
				//resclice Z so that radius measurements occur in sphere rather than elipse
				selectWindow("Loci");
				if (Isotropic == 1) {
					run("Reslice Z", "new="+W);
					closewindow("Loci");
					selectWindow("Resliced");
					rename("Loci");
				}
				
				getDimensions(dd, dd, dd, newslices, dd);
				

				
				if (LociDetType == "3D Manual Loci Selection") {
					run("Properties...", "unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
					setBatchMode("show");
					run("Set Measurements...", "centroid stack redirect=None decimal=0");
					//Stack.setChannel(LociCh);
					selectWindow("Loci");
					setTool("multipoint");
					title = "WaitForUser";
					msg = "Please select loci then click \"OK\".";
					waitForUser(title, msg);
					run("Measure");
					setBatchMode("hide");
					Xlist = newArray(0);
					Ylist = newArray(0);
					Zlist = newArray(0);
					
								
					// set counter for filtered points at 0, if point within boundary increase +1
					points_f = 0;
					for (points_i = 0; points_i < nResults; points_i++) {
						//filter the points at the beginning
				
						if ((getResult("X", points_i)) >= XY_edge && (getResult("X", points_i)) <= width-XY_edge &&
						(getResult("Y", points_i)) >= XY_edge && (getResult("Y", points_i)) <= height-XY_edge &&
						(getResult("Slice", points_i)) >= Z_edge  && (getResult("Slice", points_i)) <= newslices-Z_edge) {
								
							Xlist = append(Xlist, getResult("X", points_i));
							Ylist = append(Ylist, getResult("Y", points_i));
							Zlist = append(Zlist, getResult("Slice", points_i));
							//Ylist[points_f] = getResult("Y", points_i);
							//Zlist[points_f] = getResult("Slice", points_i);

							//print(Zlist[points_f]);
							points_f = points_f+1;
						
					
					}
					
				}
				
				print("Total points selected: "+nResults+". Total points within boundaries "+points_f+".");
				selectWindow("Loci");
				getDimensions(width, height, ChNum, slices, frames);
						
				// extract volume around centroid, detect maxima, calculate difference and update point location
													
				for (points_i = 0; points_i < points_f; points_i++) {
					// filter no longer needed as performed at start
					//if ((Xlist[points_i]) >= XY_edge && (Xlist[points_i]) <= width-XY_edge &&
					//(Ylist[points_i]) >= XY_edge && (Ylist[points_i]) <= height-XY_edge &&
					//(Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
					print(Xlist[points_i]);
					print(Ylist[points_i]);
					print(Zlist[points_i]);
					
					setSlice(Zlist[points_i]);

					// Extract volume around ROI
					leftcorner = parseInt((Xlist[points_i])-(CentroidSearchSize/2));
					topcorner = parseInt((Ylist[points_i])-(CentroidSearchSize/2));
				
					makeRectangle(leftcorner, topcorner, CentroidSearchSize, CentroidSearchSize);
					run("Duplicate...", "duplicate range="+((Zlist[points_i])-(CentroidSearchSizeZ/2))+"-"+((Zlist[points_i])+(CentroidSearchSizeZ/2))+" title=temp_region duplicate");
					run("3D Maxima Finder", "radiusxy=5 radiusz=5 noise=0");
					closewindow("temp_region");
					closewindow("peaks");
					
					Xlist[points_i] = Xlist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("X", 0)); 
					Ylist[points_i] = Ylist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("Y", 0));  
					Zlist[points_i] = Zlist[points_i] - (parseInt(CentroidSearchSizeZ/2) - getResult("Z", 0)); 
					
					run("Clear Results");

					
				//}
				}
				//create validation image
				newImage("Points", "16-bit black", width, height, slices);
				for (points_i = 0; points_i < points_f; points_i++) {
					CreateSpot(Xlist[points_i], Ylist[points_i], Zlist[points_i], 65535);
				}
				run("Merge Channels...", "c1=Points c2=Loci create keep");
				//saveAs("Tiff", LociValOut + imagename +"loci_validation.tif");
				saveAs("Tiff", ChOut + imagename +"loci_validation.tif");
				closewindow(imagename +"loci_validation.tif");
				closewindow("Points");
			
				}




				if (LociDetType == "2D Manual Loci Selection") {
					run("Properties...", "unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
					run("Z Project...", "projection=[Max Intensity]");
					selectWindow("MAX_Loci");
					setBatchMode("show");
					run("Set Measurements...", "centroid stack redirect=None decimal=0");
					//Stack.setChannel(LociCh);
					
					setTool("multipoint");
					title = "WaitForUser";
					msg = "Please select loci then click \"OK\".";
					waitForUser(title, msg);
					run("Measure");
					setBatchMode("hide");
					closewindow("MAX_Loci");
					Xlist = newArray(0);
					Ylist = newArray(0);
					Zlist = newArray(0);
					// set counter for filtered points at 0, if point within boundary increase +1
					points_f = 0;
					for (points_i = 0; points_i < nResults; points_i++) {
						//filter the points at the beginning
											
						if ((getResult("X", points_i)) >= XY_edge && (getResult("X", points_i)) <= width-XY_edge &&
						(getResult("Y", points_i)) >= XY_edge && (getResult("Y", points_i)) <= height-XY_edge)	 {
						
						
							Xlist = append(Xlist, getResult("X", points_i));
							Ylist = append(Ylist, getResult("Y", points_i));
							Zlist = append(Zlist, 5);
							//Ylist[points_f] = getResult("Y", points_i);
							//Zlist[points_f] = getResult("Slice", points_i);

							//print(Zlist[points_f]);
							points_f = points_f+1;
						
					
					}
					
				}
				print("Total points selected: "+nResults+". Total points within boundaries "+points_f+".");
				
				selectWindow("Loci");
				getDimensions(width, height, ChNum, slices, frames);
						
					// extract volume around centroid, detect maxima, calculate difference and update point location
													
					for (points_i = 0; points_i < points_f; points_i++) {
						// filter no longer needed as performed at start
						//if ((Xlist[points_i]) >= XY_edge && (Xlist[points_i]) <= width-XY_edge &&
						//(Ylist[points_i]) >= XY_edge && (Ylist[points_i]) <= height-XY_edge &&
						//(Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
						
						// Extract volume around ROI
						leftcorner = parseInt((Xlist[points_i])-(CentroidSearchSize/2));
						topcorner = parseInt((Ylist[points_i])-(CentroidSearchSize/2));
					
						makeRectangle(leftcorner, topcorner, CentroidSearchSize, CentroidSearchSize);
						run("Duplicate...", "duplicate title=temp_region duplicate");
						run("3D Maxima Finder", "radiusxy=5 radiusz=5 noise=0");
						closewindow("temp_region");
						closewindow("peaks");
						
						Xlist[points_i] = Xlist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("X", 0)); 
						Ylist[points_i] = Ylist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("Y", 0));  
						Zlist[points_i] = getResult("Z", 0); 
						run("Clear Results");

						print(Zlist[points_i]);
					//}
					}
					//create validation image
					newImage("Points", "16-bit black", width, height, slices);
					for (points_i = 0; points_i < points_f; points_i++) {
						if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
							CreateSpot(Xlist[points_i], Ylist[points_i], Zlist[points_i], 65535);
						}
					}
					run("Merge Channels...", "c1=Points c2=Loci create keep");
					//saveAs("Tiff", LociValOut + imagename +"loci_validation.tif");
					saveAs("Tiff", ChOut + imagename +"loci_validation.tif");
					closewindow(imagename +"loci_validation.tif");
					closewindow("Points");

				}

   
				
				// Process First Cell Channel

				// Run Pre Filter - note customised for # of CPUs - actually not necessary. Binarize Loci channel then multiply with Raw data to get clean loci
	
				//run("3D Fast Filters","filter=TopHat radius_x_pix=2.0 radius_y_pix=2.0 radius_z_pix=2.0 Nb_cpus=12");
				//selectWindow("3D_TopHat");

				// Enhance and clean

				//run("Unsharp Mask...", "radius=5 mask=0.70 stack");

				//RemoveOutliersFilter("Loci");
				
				//run("Median 3D...", "x=5 y=5 z=5");
				if (LociDetType == "Automatic Detection") {
					run("Duplicate...", "title=Loci-Mask duplicate");
					selectWindow("Loci-Mask");
					run("Properties...", "unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
							
					setAutoThreshold("RenyiEntropy dark no-reset stack"); // before filtering and unsharp were used
					//run("Auto Threshold", "method=Li white stack use_stack_histogram");
					run("Convert to Mask", "method=RenyiEntropy background=Dark black");
					
					run("3D OC Options", "centre_of_mass dots_size=5 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
					run("3D Objects Counter", "threshold=128 slice=30 min.=2 max.=6877872 statistics");
					Table.rename("Statistics for Loci-Mask", "Results");
					closewindow("Loci-Mask");

					Xlist = newArray(0);
					Ylist = newArray(0);
					Zlist = newArray(0);
					points_f = 0;
					for (points_i = 0; points_i < nResults; points_i++) {
						//filter the points at the beginning
												
						if ((getResult("XM", points_i)) >= XY_edge && (getResult("XM", points_i)) <= width-XY_edge &&
						(getResult("YM", points_i)) >= XY_edge && (getResult("YM", points_i)) <= height-XY_edge &&
						(getResult("ZM", points_i)) >= Z_edge  && (getResult("ZM", points_i)) <= newslices-Z_edge) {

							Xlist = append(Xlist, parseInt(getResult("XM", points_i)));
							Ylist = append(Ylist, parseInt(getResult("YM", points_i)));
							Zlist = append(Zlist, parseInt(getResult("ZM", points_i)));
							//Ylist[points_f] = getResult("Y", points_i);
							//Zlist[points_f] = getResult("Slice", points_i);

							//print(Zlist[points_f]);
							points_f = points_f+1;
						
					
					}
					}
				
					print("Total points selected: "+nResults+". Total points within boundaries "+points_f+".");

					selectWindow("Loci");
					getDimensions(width, height, ChNum, slices, frames);
					
					for (points_i = 0; points_i < points_f; points_i++) {
						// filter no longer needed as performed at start
						//if ((Xlist[points_i]) >= XY_edge && (Xlist[points_i]) <= width-XY_edge &&
						//(Ylist[points_i]) >= XY_edge && (Ylist[points_i]) <= height-XY_edge &&
						//(Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
					
						
						setSlice(Zlist[points_i]);

						// Extract volume around ROI
						leftcorner = parseInt((Xlist[points_i])-(CentroidSearchSize/2));
						topcorner = parseInt((Ylist[points_i])-(CentroidSearchSize/2));
					
						makeRectangle(leftcorner, topcorner, CentroidSearchSize, CentroidSearchSize);
						run("Duplicate...", "duplicate range="+((Zlist[points_i])-(CentroidSearchSizeZ/2))+"-"+((Zlist[points_i])+(CentroidSearchSizeZ/2))+" title=temp_region duplicate");
						run("3D Maxima Finder", "radiusxy=5 radiusz=5 noise=0");
						closewindow("temp_region");
						closewindow("peaks");
						
						Xlist[points_i] = Xlist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("X", 0)); 
						Ylist[points_i] = Ylist[points_i] - (parseInt(CentroidSearchSize/2) - getResult("Y", 0));  
						Zlist[points_i] = Zlist[points_i] - (parseInt(CentroidSearchSizeZ/2) - getResult("Z", 0)); 
						
						run("Clear Results");

						print(Zlist[points_i]);
					//}
					}
					//create validation image
					newImage("Points", "16-bit black", width, height, newslices);
					for (points_i = 0; points_i < points_f; points_i++) {
						CreateSpot(Xlist[points_i], Ylist[points_i], Zlist[points_i], 65535);
					}
					run("Merge Channels...", "c1=Points c2=Loci create keep");
					//saveAs("Tiff", LociValOut + imagename +"loci_validation.tif");
					saveAs("Tiff", ChOut + imagename +"loci_validation.tif");
					closewindow(imagename +"loci_validation.tif");
					closewindow("Points");
		
				}
				

				
error_count = 0;
error_countv = 0;

// Measure each point in DistCh1
point_counter = 0;


if (Xlist.length > 0) {

				

//---Process Channel of Interest---Measuring from real Loci locations--------------
			
				//resclice Z so that radius measurements occur in sphere rather than elipse
				selectImage(RDCh1ID);
				rename("RDCh1");
				if (Isotropic == 1) {
					run("Reslice Z", "new="+W);
					closewindow("RDCh1");
					selectWindow("Resliced");
					rename("RDCh1");
				}
				
				//Table for Measured Intensities
				MeasuredIntTable1 = "MeasuredInt_of_Radial_Distributions"; 
				MeasuredIntTable2 = "["+MeasuredIntTable1+"]"; 
				MeasuredIntTableTitle=MeasuredIntTable2; 
				run("New... ", "name="+MeasuredIntTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?			

				//Create Table Headings
				TableHeadings = "\\Headings:Loci";
				for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
					TableHeadings = TableHeadings + "\t"+array_i;
				}
				
				print(MeasuredIntTableTitle,TableHeadings); 

				
				//Table for Normalized Intensities			
				NormIntTable1 = "Summary_of_Radial_Distributions"; 
				NormIntTable2 = "["+NormIntTable1+"]"; 
				NormIntTableTitle=NormIntTable2; 
				run("New... ", "name="+NormIntTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?
				//print(NormIntTableTitle,"\\Headings:Loci\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46"); 
				print(NormIntTableTitle, TableHeadings); 


				//Table for Loci Summary			
				SummaryTable1 = "Summary_of_Loci"; 
				SummaryTable2 = "["+SummaryTable1+"]"; 
				SummaryTableTitle=SummaryTable2; 
				run("New... ", "name="+SummaryTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?
				print(SummaryTableTitle,"\\Headings:Loci\tX\tY\tZ\tCenter_Intensity\t250um_Mean_Intensity\tNucleus_Mean_Intensity\tRD_Min_Intensity\tNormalized_Center_Intensity\tCenter_to_Nucleus_Ratio\t250um_Sphere_to_Nucleus_Ratio"); 
			


				
				//create an empty array for normalized intensity at centre
				Intlist = newArray(points_f);
				//create an array for the mean intensity of the nucleus
				MeanIntNucleusList = newArray(points_f);
				//create an array for the mean intensity at 250um divided by the nucleus mean itensity * 100 for plotting ratio
				MeanInt250umRatio = newArray(points_f);

				
				//process each point
				run("Clear Results");
				closewindow("Results");
				for (points_i = 0; points_i < points_f; points_i++) {
					if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
					point_counter ++;
					//make sure points are not within the boundary
					// filter no longer needed as performed at start
					//if ((Xlist[points_i]) >= XY_edge && (Xlist[points_i]) <= width-XY_edge &&
					//	(Ylist[points_i]) >= XY_edge && (Ylist[points_i]) <= height-XY_edge &&
					//	(Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
							
						
					selectWindow("RDCh1");
					setSlice(Zlist[points_i]);

					// Extract ROI and save to combine later (should make this a function) MakeBoxROI(x,y,xlen,ylen)
				
					//print(testing);
					leftcorner = parseInt((Xlist[points_i])-(ExtractROIsize/2));
					if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
						leftcorner = 0;
					}
					topcorner = parseInt((Ylist[points_i])-(ExtractROIsize/2));
					if (topcorner < 0) {
						//boxheight = boxheight +(topcorner*-1);
						topcorner = 0;
					}
					makeRectangle(leftcorner, topcorner, ExtractROIsize, ExtractROIsize);
					
					run("Duplicate...", " ");
					run("Canvas Size...", "width="+ExtractROIsize+" height="+ExtractROIsize+" position=Center zero");
					ROIfilename = points_i + 10000;
					saveAs("Tiff", MeasuredOut + ROIfilename +".tif");
					close();


					
					//make a square to measure the mean intensity of the nucleus 
					run("Set Measurements...", "mean limit redirect=None decimal=0");
					
					leftcorner = parseInt((Xlist[points_i])-(box4um/2));
					if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
						leftcorner = 0;
					}
					topcorner = parseInt((Ylist[points_i])-(box4um/2));
					if (topcorner < 0) {
						//boxheight = boxheight +(topcorner*-1);
						topcorner = 0;
					}
					makeRectangle(leftcorner, topcorner, box4um, box4um);
					setAutoThreshold("Otsu dark no-reset");
					run("Measure");
					MeanIntNucleusList[points_i] = getResult("Mean", 0);
					

					run("Clear Results");
		
					//place point and measure radial dist
					makePoint(Xlist[points_i], Ylist[points_i]);	
					
					run("3D Radial Distribution", "radius_max="+ RDRadSteps +" measure=Mean");
				
					if (isOpen("Radial distribution") == 1){
						Plot.getValues(xpoints, ypoints);
						normval = ypoints[RDRadSteps];
					} else {
						error_count = error_count+1;
						closewindow("Exception");
						xpoints = newArray(RDRadSteps*2+1);
						ypoints = newArray(RDRadSteps*2+1);
						normval = 0;
						//closewindow("Exception");
						//print("Reattempt RD+10 analysis");
						//run("3D Radial Distribution", "radius_max="+ RDRadSteps+10 +" measure=Mean");
						//if (isOpen("Radial distribution") == 1){
						//	Plot.getValues(xpoints, ypoints);
						//	normval = ypoints[RDRadSteps+10];
						//} else {
						//	closewindow("Exception");
						//	print("creating crop region for RD analysis");
							// Extract volume around ROI
						//	leftcorner = parseInt((Xlist[points_i])-(RadExtract/2));
						//	topcorner = parseInt((Ylist[points_i])-(RadExtract/2));
						
						//	makeRectangle(leftcorner, topcorner, RadExtract, RadExtract);
						//	run("Duplicate...", "duplicate range="+((Zlist[points_i])-(RadExtract/2))+"-"+((Zlist[points_i])+(RadExtract/2))+" title=temp_region duplicate");	
						//	selectWindow("temp_region");
						//	setSlice(RadExtract/2);
						//	makePoint(RadExtract/2, RadExtract/2);					
						//	run("3D Radial Distribution", "radius_max="+ RDRadSteps+10 +" measure=Mean");
						//	if (isOpen("Radial distribution") == 1){
						//		Plot.getValues(xpoints, ypoints);
						//		normval = ypoints[RDRadSteps+10];
						//	}
						//closewindow("Exception");
						//closewindow("temp_region");
						
					//}
					}
					
				
					
					sphere250umInt = 0;
					// measure sphere250umInt
					for (spherepoints = 0; spherepoints < sphere250; spherepoints++){
						sphere250umInt = sphere250umInt + ypoints[RDRadSteps + spherepoints];
					}
					sphere250umInt = sphere250umInt/sphere250;
					
					tab = "\t";
					//normalise values to 1
					NormDistArray = newArray(RDRadSteps+1);
					// Create Normalized array
					for (plot_i=RDRadSteps; plot_i<((RDRadSteps*2)+1); plot_i++) {
						NormDistArray[plot_i-RDRadSteps] = ypoints[plot_i]/normval;
					}

					// Array with background removed
					Array.getStatistics(NormDistArray, NDmin, dummy, dummy, dummy);
						for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						NormDistArray[array_i] = NormDistArray[array_i] - NDmin;
					}

					Intlist[points_i] = parseInt((NormDistArray[0]*100));
		
					
					// generate measured intensities table line
					printlineint = toString(points_i+1);
					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printlineint = printlineint + tab + ypoints[array_i+RDRadSteps];
					}
					print(MeasuredIntTableTitle, printlineint);


					
					// generate profile table line
					printline = toString(points_i+1);

					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printline = printline + tab + NormDistArray[array_i];
					}
					//for (plot_i=46; plot_i<93; plot_i++) {
						//	printline = printline + tab + ((ypoints[plot_i])/normval);
					//}
					closewindow("Radial distribution");
					print(NormIntTableTitle, printline);

					MeanInt250umRatio[points_i] = (sphere250umInt/MeanIntNucleusList[points_i]*100);

					// Create summary table line
				
					printlineint = toString(points_i+1);
					printlineint = printlineint + tab + Xlist[points_i] + tab + Ylist[points_i] + tab + Zlist[points_i] + tab + normval + tab + sphere250umInt + tab + MeanIntNucleusList[points_i] + tab + (normval*NDmin) + tab + NormDistArray[0] + tab + (normval/MeanIntNucleusList[points_i]) + tab + (sphere250umInt/MeanIntNucleusList[points_i]); 
					print(SummaryTableTitle, printlineint);

				//}
				}
				}

			
				//save table raw intensities
				selectWindow(MeasuredIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Raw_Radial_Distribution_COI_Ch_"+DistCh1+".csv]");
				closewindow(MeasuredIntTableTitle);				
				
				
				//save table Normalized
						
				selectWindow(NormIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Normalized_Radial_Distribution_COI_Ch_"+DistCh1+".csv]");
				closewindow(NormIntTableTitle);

				//save summary table
						
				selectWindow(SummaryTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Summary_of_Loci_Measurements_COI_Ch_"+DistCh1+".csv]");
				closewindow(SummaryTableTitle);
				
				// generate new image contain colored loci based on normalised int
				newImage("Coded_Loci", "16-bit black", width, height, newslices);
				
				for (points_i = 0; points_i < points_f; points_i++) {
					if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
						selectWindow("Coded_Loci");
						CreateSphereAdditive3px(Xlist[points_i], Ylist[points_i], Zlist[points_i], Intlist[points_i]);
					}
				}
				
				run("Rainbow RGB");
				saveAs("Tiff", ChOut + imagename +"_Loci_Intensity_COI_Ch_"+DistCh1+".tif");
				close();

				
				// generate new image contain colored loci based on ratio of sphere to nucleus
				
				newImage("Coded_Loci", "16-bit black", width, height, newslices);
				
				for (points_i = 0; points_i < points_f; points_i++) {
					if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
						selectWindow("Coded_Loci");
						CreateSphereAdditive3px(Xlist[points_i], Ylist[points_i], Zlist[points_i], MeanInt250umRatio[points_i]);
					}
				}
				
				run("Rainbow RGB");
				saveAs("Tiff", ChOut + imagename +"_Loci_250_Sphere_Nucleus_Ratio__COI_Ch_"+DistCh1+".tif");
				close();



				
//------Process Channel of Interest--- MEASURING RANDOMLY SEEDED FAUX LOCI--------------
			
				selectWindow("RDCh1");

				//Table for Measured Intensities
				MeasuredIntTable1 = "MeasuredInt_of_Radial_Distributions"; 
				MeasuredIntTable2 = "["+MeasuredIntTable1+"]"; 
				MeasuredIntTableTitle=MeasuredIntTable2; 
				run("New... ", "name="+MeasuredIntTable2+" type=Table"); 
				
				//Create Table Headings
				TableHeadings = "\\Headings:Loci";
				for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
					TableHeadings = TableHeadings + "\t"+array_i;
				}
				
				print(MeasuredIntTableTitle,TableHeadings); 

				//Table for Normalized Intensities			
				NormIntTable1 = "Summary_of_Radial_Distributions"; 
				NormIntTable2 = "["+NormIntTable1+"]"; 
				NormIntTableTitle=NormIntTable2; 
				run("New... ", "name="+NormIntTable2+" type=Table"); 
				print(NormIntTableTitle, TableHeadings); 

				//Table for Loci Summary			
				SummaryTable1 = "Summary_of_Loci"; 
				SummaryTable2 = "["+SummaryTable1+"]"; 
				SummaryTableTitle=SummaryTable2; 
				run("New... ", "name="+SummaryTable2+" type=Table"); 
				print(SummaryTableTitle,"\\Headings:Loci\tX\tY\tZ\tCenter_Intensity\t250um_Mean_Intensity\tNucleus_Mean_Intensity\tRD_Min_Intensity\tNormalized_Center_Intensity\tCenter_to_Nucleus_Ratio\t250um_Sphere_to_Nucleus_Ratio"); 
			
				//create an empty array for normalized intensity at centre
				Intlist = newArray(points_f);
				//create an array for the mean intensity of the nucleus
				MeanIntNucleusList = newArray(points_f);
				//create an array for the mean intensity at 250um divided by the nucleus mean itensity * 100 for plotting ratio
				MeanInt250umRatio = newArray(points_f);
				
				
				run("Clear Results");
				closewindow("Results");

				//CREATE AND ANALYZE A RANDOM LOCI with background above Minimum specified				
				for (points_i = 0; points_i < point_counter; points_i++) {
					selectWindow("RDCh1");
					for (attempts=1; attempts<100; attempts++) {
						run("Set Measurements...", "mean limit redirect=None decimal=0");
						RanX = parseInt(random() * (width-(XY_edge*2))) + XY_edge;
						RanY = parseInt(random() * (height-(XY_edge*2)))+ XY_edge;
						RanZ = parseInt(random() * (slices-(Z_edge*2))) + Z_edge;

						setSlice(RanZ);
						makePoint(RanX, RanY);
						run("Measure");
						if (getResult("Mean", 0) >= SeedMinInt) {
							attempts = 100;
						}
						run("Clear Results");
					}
							


				
					// Extract ROI and save to combine later (should make this a function) MakeBoxROI(x,y,xlen,ylen)
				
					//leftcorner = parseInt((RanX)-(ExtractROIsize/2));
					//if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
					//	leftcorner = 0;
					//}
					//topcorner = parseInt((RanY)-(ExtractROIsize/2));
					//if (topcorner < 0) {
					//	//boxheight = boxheight +(topcorner*-1);
					//	topcorner = 0;
					//}
					//makeRectangle(leftcorner, topcorner, ExtractROIsize, ExtractROIsize);
					
					//run("Duplicate...", " ");
					//run("Canvas Size...", "width="+ExtractROIsize+" height="+ExtractROIsize+" position=Center zero");
					//ROIfilename = points_i + 10000;
					//saveAs("Tiff", MeasuredOut + ROIfilename +".tif");
					//close();

					//make a square to measure the mean intensity of the nucleus 
					//run("Set Measurements...", "mean limit redirect=None decimal=0");
					
					leftcorner = parseInt((RanX)-(box4um/2));
					if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
						leftcorner = 0;
					}
					topcorner = parseInt((RanY)-(box4um/2));
					if (topcorner < 0) {
						//boxheight = boxheight +(topcorner*-1);
						topcorner = 0;
					}
					makeRectangle(leftcorner, topcorner, box4um, box4um);
					setAutoThreshold("Otsu dark no-reset");
					run("Measure");
					MeanIntNucleusList[points_i] = getResult("Mean", 0);

					run("Clear Results");
		
					//place point and measure radial dist
					makePoint(RanX, RanY);	
					
					run("3D Radial Distribution", "radius_max="+ RDRadSteps +" measure=Mean");
				
					if (isOpen("Radial distribution") == 1){
						Plot.getValues(xpoints, ypoints);
						normval = ypoints[RDRadSteps];
					} else {
						error_count = error_count+1;
						closewindow("Exception");
						xpoints = newArray(RDRadSteps*2+1);
						ypoints = newArray(RDRadSteps*2+1);
						normval = 0;
			
					}
					
				
					
					sphere250umInt = 0;
					// measure sphere250umInt
					for (spherepoints = 0; spherepoints < sphere250; spherepoints++){
						sphere250umInt = sphere250umInt + ypoints[RDRadSteps + spherepoints];
					}
					sphere250umInt = sphere250umInt/sphere250;
					
					tab = "\t";
					//normalise values to 1
					NormDistArray = newArray(RDRadSteps+1);
					// Create Normalized array
					for (plot_i=RDRadSteps; plot_i<((RDRadSteps*2)+1); plot_i++) {
						NormDistArray[plot_i-RDRadSteps] = ypoints[plot_i]/normval;
					}

					// Array with background removed
					Array.getStatistics(NormDistArray, NDmin, dummy, dummy, dummy);
						for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						NormDistArray[array_i] = NormDistArray[array_i] - NDmin;
					}

					Intlist[points_i] = parseInt((NormDistArray[0]*100));
		
					
					// generate measured intensities table line
					printlineint = toString(points_i+1);
					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printlineint = printlineint + tab + ypoints[array_i+RDRadSteps];
					}
					print(MeasuredIntTableTitle, printlineint);


					
					// generate profile table line
					printline = toString(points_i+1);

					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printline = printline + tab + NormDistArray[array_i];
					}
					//for (plot_i=46; plot_i<93; plot_i++) {
						//	printline = printline + tab + ((ypoints[plot_i])/normval);
					//}
					closewindow("Radial distribution");
					print(NormIntTableTitle, printline);

					MeanInt250umRatio[points_i] = (sphere250umInt/MeanIntNucleusList[points_i]*100);

					// Create summary table line
				
					printlineint = toString(points_i+1);
					printlineint = printlineint + tab + RanX + tab + RanY + tab + RanZ + tab + normval + tab + sphere250umInt + tab + MeanIntNucleusList[points_i] + tab + (normval*NDmin) + tab + NormDistArray[0] + tab + (normval/MeanIntNucleusList[points_i]) + tab + (sphere250umInt/MeanIntNucleusList[points_i]); 
					print(SummaryTableTitle, printlineint);

				//}
				
				}

			
				//save table raw intensities
				selectWindow(MeasuredIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Raw_Radial_Distribution_COI_FauxLoci_Ch_"+DistCh1+".csv]");
				closewindow(MeasuredIntTableTitle);				
				
				
				//save table Normalized
						
				selectWindow(NormIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Normalized_Radial_Distribution_COI_FauxLoci_Ch_"+DistCh1+".csv]");
				closewindow(NormIntTableTitle);

				//save summary table
						
				selectWindow(SummaryTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Summary_of_Loci_Measurements_COI_FauxLoci_Ch_"+DistCh1+".csv]");
				closewindow(SummaryTableTitle);
				
				// generate new image contain colored loci based on normalised int
				//newImage("Coded_Loci", "16-bit black", width, height, newslices);
				//
				//for (points_i = 0; points_i < points_f; points_i++) {
				//	if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
				//		selectWindow("Coded_Loci");
				//		CreateSphereAdditive3px(Xlist[points_i], Ylist[points_i], Zlist[points_i], Intlist[points_i]);
				//	}
				//}
				
				//run("Rainbow RGB");
				//saveAs("Tiff", ChOut + imagename +"_Loci_intensity_COI_FauxLoci_Ch_"+DistCh1+".tif");
				//close();

				
				// generate new image contain colored loci based on ratio of sphere to nucleus
				
				//newImage("Coded_Loci", "16-bit black", width, height, newslices);
				
				//for (points_i = 0; points_i < points_f; points_i++) {
				//	if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
				//		selectWindow("Coded_Loci");
				//		CreateSphereAdditive3px(Xlist[points_i], Ylist[points_i], Zlist[points_i], MeanInt250umRatio[points_i]);
				//	}
				//}
				
				//run("Rainbow RGB");
				//saveAs("Tiff", ChOut + imagename +"_loci_250_Sphere_Nucleus_Ratio_COI_FauxLoci_Ch_"+DistCh1+".tif");
				//close();

								


	


//----ADDED FOR VALIDATION ---- Process Loci channel in same way

//cleaning up images and names, remove later
selectWindow("Loci");		
rename("Loci-Raw");	
				
				selectWindow("Loci-Raw");
				//Table for Measured Intensities
				MeasuredIntTable1 = "MeasuredInt_of_Radial_Distributions"; 
				MeasuredIntTable2 = "["+MeasuredIntTable1+"]"; 
				MeasuredIntTableTitle=MeasuredIntTable2; 
				run("New... ", "name="+MeasuredIntTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?
				//print(MeasuredIntTableTitle,"\\Headings:Loci\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46"); 
				print(MeasuredIntTableTitle, TableHeadings); 
				
	
				//Table for Normalized Intensities			
				NormIntTable1 = "Summary_of_Radial_Distributions"; 
				NormIntTable2 = "["+NormIntTable1+"]"; 
				NormIntTableTitle=NormIntTable2; 
				run("New... ", "name="+NormIntTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?
				//print(NormIntTableTitle,"\\Headings:Loci\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46"); 
				print(NormIntTableTitle, TableHeadings); 

				
				//Table for Loci Summary			
				SummaryTable1 = "Summary_of_Loci"; 
				SummaryTable2 = "["+SummaryTable1+"]"; 
				SummaryTableTitle=SummaryTable2; 
				run("New... ", "name="+SummaryTable2+" type=Table"); 
				//Titles only going to 46 (original radius) - if wanting to plot longer than need to make a loop that keeps adding slashed seperated numbers. is this necessary?
				print(SummaryTableTitle,"\\Headings:Loci\tX\tY\tZ\tCenter_Intensity\t250um_Mean_Intensity\tNucleus_Mean_Intensity\tRD_Min_Intensity\tNormalized_Center_Intensity\tCenter_to_Nucleus_Ratio\t250um_Sphere_to_Nucleus_Ratio"); 
			


				
				//create an empty array for normalized intensity at centre
				Intlist = newArray(points_f);
				//create an array for the mean intensity of the nucleus
				MeanIntNucleusList = newArray(points_f);
				//create an array for the mean intensity at 250um divided by the nucleus mean itensity * 100 for plotting ratio
				MeanInt250umRatio = newArray(points_f);

				run("Clear Results");
				closewindow("Results");
				//process each point
				for (points_i = 0; points_i < points_f; points_i++) {
					if ((Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
					// filter no longer needed as performed at start
					//if ((Xlist[points_i]) >= XY_edge && (Xlist[points_i]) <= width-XY_edge &&
					//	(Ylist[points_i]) >= XY_edge && (Ylist[points_i]) <= height-XY_edge &&
					//	(Zlist[points_i]) >= Z_edge  && (Zlist[points_i]) <= slices-Z_edge) {
							
						
					selectWindow("Loci-Raw");
					setSlice(Zlist[points_i]);

					// Extract ROI and save to combine later (should make this a function) MakeBoxROI(x,y,xlen,ylen)
					leftcorner = parseInt((Xlist[points_i])-(ExtractROIsize/2));
					if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
						leftcorner = 0;
					}
					topcorner = parseInt((Ylist[points_i])-(ExtractROIsize/2));
					if (topcorner < 0) {
						//boxheight = boxheight +(topcorner*-1);
						topcorner = 0;
					}
					makeRectangle(leftcorner, topcorner, ExtractROIsize, ExtractROIsize);
					run("Duplicate...", " ");
					run("Canvas Size...", "width="+ExtractROIsize+" height="+ExtractROIsize+" position=Center zero");
					ROIfilename = points_i + 10000;
					saveAs("Tiff", LociOut + ROIfilename +".tif");
					close();
					
					//make a squre to measure the mean intensity of the nucleus 
					run("Set Measurements...", "mean limit redirect=None decimal=0");
					
					leftcorner = parseInt((Xlist[points_i])-(box4um/2));
					if (leftcorner < 0) {
						//boxwidth = boxwidth +(leftcorner*-1);
						leftcorner = 0;
					}
					topcorner = parseInt((Ylist[points_i])-(box4um/2));
					if (topcorner < 0) {
						//boxheight = boxheight +(topcorner*-1);
						topcorner = 0;
					}
					makeRectangle(leftcorner, topcorner, box4um, box4um);
					setAutoThreshold("Otsu dark no-reset");
					run("Measure");
					MeanIntNucleusList[points_i] = getResult("Mean", 0);
					run("Clear Results");

					
					//place point and measure radial dist
					makePoint(Xlist[points_i], Ylist[points_i]);					
					run("3D Radial Distribution", "radius_max="+ RDRadSteps +" measure=Mean");
					
					if (isOpen("Radial distribution") == 1){
						Plot.getValues(xpoints, ypoints);
						normval = ypoints[RDRadSteps];
					} else {
						error_countv = error_count+1;
						closewindow("Exception");
						xpoints = newArray(RDRadSteps*2+1);
						ypoints = newArray(RDRadSteps*2+1);
						normval = 0;
						
						//print("Reattempt RD+10 analysis");
						//run("3D Radial Distribution", "radius_max="+ RDRadSteps+10 +" measure=Mean");
						//if (isOpen("Radial distribution") == 1){
						//	Plot.getValues(xpoints, ypoints);
						//	normval = ypoints[RDRadSteps+10];
						//} else {
						//	closewindow("Exception");
						//	print("creating crop region for RD analysis");
							// Extract volume around ROI
						//	leftcorner = parseInt((Xlist[points_i])-(RadExtract/2));
						//	topcorner = parseInt((Ylist[points_i])-(RadExtract/2));
						
						//	makeRectangle(leftcorner, topcorner, RadExtract, RadExtract);
						//	run("Duplicate...", "duplicate range="+((Zlist[points_i])-(RadExtract/2))+"-"+((Zlist[points_i])+(RadExtract/2))+" title=temp_region duplicate");	
						//	selectWindow("temp_region");
						//	setSlice(RadExtract/2);
						//	makePoint(RadExtract/2, RadExtract/2);					
						//	run("3D Radial Distribution", "radius_max="+ RDRadSteps+10 +" measure=Mean");
						//	if (isOpen("Radial distribution") == 1){
						//		Plot.getValues(xpoints, ypoints);
						//		normval = ypoints[RDRadSteps+10];
						//	}
						//closewindow("Exception");
						//closewindow("temp_region");
						
					//}
					}
					
					
					
					sphere250umInt = 0;
					// measure sphere250umInt
					for (spherepoints = 0; spherepoints < sphere250; spherepoints++){
						sphere250umInt = sphere250umInt + ypoints[RDRadSteps + spherepoints];
					}
					sphere250umInt = sphere250umInt/sphere250;
					
					tab = "\t";
					//normalise values to 1
					NormDistArray = newArray(RDRadSteps+1);
					// Create Normalized array
					for (plot_i=RDRadSteps; plot_i<((RDRadSteps*2)+1); plot_i++) {
						NormDistArray[plot_i-RDRadSteps] = ypoints[plot_i]/normval;
					}

					// Array with background removed
					Array.getStatistics(NormDistArray, NDmin, dummy, dummy, dummy);
						for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						NormDistArray[array_i] = NormDistArray[array_i] - NDmin;
					}

					Intlist[points_i] = parseInt((NormDistArray[0]*100));
		
					
					// generate measured intensities table line
					printlineint = toString(points_i+1);
					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printlineint = printlineint + tab + ypoints[array_i+RDRadSteps];
					}
					print(MeasuredIntTableTitle, printlineint);


					
					// generate profile table line
					printline = toString(points_i+1);

					for (array_i=0; array_i<(RDRadSteps+1); array_i++) {
						printline = printline + tab + NormDistArray[array_i];
					}
					//for (plot_i=46; plot_i<93; plot_i++) {
						//	printline = printline + tab + ((ypoints[plot_i])/normval);
					//}
					closewindow("Radial distribution");
					print(NormIntTableTitle, printline);

					MeanInt250umRatio[points_i] = (sphere250umInt/MeanIntNucleusList[points_i]*100);

					// Create summary table line
					printlineint = toString(points_i+1);
					printlineint = printlineint + tab + Xlist[points_i] + tab +Ylist[points_i] + tab + Zlist[points_i] + tab + normval + tab + sphere250umInt + tab + MeanIntNucleusList[points_i] + tab +(normval*NDmin) + tab + NormDistArray[0] + tab + (normval/MeanIntNucleusList[points_i]) + tab + (sphere250umInt/MeanIntNucleusList[points_i]); 
					print(SummaryTableTitle, printlineint);

					
	
				//}
				}
				}

				
				//save table raw intensities
				selectWindow(MeasuredIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Raw_Radial_Distribution_Loci_Ch_"+LociCh+".csv]");
				closewindow(imagename  + "_Raw_Radial_Distribution_measurements_Ch_"+LociCh+".csv");				
				
				
				//save table Normalized
						
				selectWindow(NormIntTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Normalized_Radial_Distribution_Loci_Ch_"+LociCh+".csv]");
				closewindow(imagename  + "_Normalized_Radial_Distribution_measurements_Ch_"+LociCh+".csv");

				//save summary table
						
				selectWindow(SummaryTable1);
				run("Text...", "save=["+ ChOut + imagename  + "_Summary_of_Loci_Measurements_LociCh_Ch_"+LociCh+".csv]");
				closewindow(imagename  + "_Summary_of_Loci_Measurements_Ch_"+LociCh+".csv");

				//Create Zstack of points and montage
				run("Image Sequence...", "open=["+MeasuredOut+ "10000.tif] sort");
				rename("Measured");
				run("Image Sequence...", "open=["+LociOut+ "10000.tif] sort");
				rename("Loci");
				run("Merge Channels...", "c2=Loci c6=Measured create");
				saveAs("Tiff", ChOut + imagename +"_Loci_Locations.tif");
				
				close();
				
				DeleteDir(MeasuredOut);
				DeleteDir(LociOut);
			


          	CloseOpenWindows();
          	closewindow("Loci-Raw");
          	closewindow("C1-Raw");
          	closewindow("C3-Raw");
          	collectGarbage(20, 4);
        }
        
		DeleteDir(TempOut);
		endtime = getTime();
		dif = (endtime-starttime)/1000;
		print("  ");
		print("----------------------------------------------------------------");
		print("Loci Radial Distribution Analysis Complete. Processing time =", (dif/60), "minutes. ", (dif/files.length), "seconds per image.");
		print("----------------------------------------------------------------");
		print("3D Radial Distribution Error on Full Res Image "+error_count+" times, and validation image "+error_countv+" times. Line substitute with zeros on these occasions.");
		selectWindow("Log");
		run("Text...", "save=["+ ChOut + imagename  + "_Loci_Analysis_Log.txt]");
		closewindow("RDCh1");
		
		
		}
	
	
   }
}
}


function closewindow(windowname) {
	if (isOpen(windowname)) { 
		 selectWindow(windowname); 
 		run("Close"); 
  	} 
}
function collectGarbage(slices, itr){
	setBatchMode(false);
	wait(1000);
	for(i=0; i<itr; i++){
		wait(50*slices);
		run("Collect Garbage");
		call("java.lang.System.gc");
		}
	setBatchMode(true);
}

function RemoveOutliersFilter(imagename) {
	selectWindow(imagename);
	getDimensions(w2, h2, c2, slices, f2);
	rename("ROFimage");
	run("Duplicate...", "title=bgstack duplicate");
	run("Z Project...", "projection=[Max Intensity]");
	for(i=0; i<30; i++) {	
		run("Remove Outliers...", "radius=3 threshold=0 which=Bright");
		//run("Morphological Filters", "operation=Dilation element=Disk radius=3");
	}
	run("Morphological Filters", "operation=Dilation element=Disk radius=10");
	rename("morph");
	selectWindow("MAX_bgstack");
	close();
	selectWindow("morph");
	rename("MAX_bgstack");
	
	for (i=0; i<slices; i++) {
		selectWindow("MAX_bgstack");
		run("Select All");
		run("Copy");
		selectWindow("bgstack");
		setSlice(i+1);
		run("Paste");
	}
	selectWindow("MAX_bgstack");
	close();
	
	imageCalculator("Subtract create stack", "ROFimage","bgstack");
	selectWindow("bgstack");
	close();
	selectWindow("ROFimage");
	close();
	selectWindow("Result of ROFimage");
	rename(imagename);
}	
function OSBSFilter(imagename, radius, iterations) {
	selectWindow(imagename);
	getDimensions(w2, h2, c2, slices, f2);
	rename("ROFimage");
	run("Duplicate...", "title=bgstack duplicate");
	run("Z Project...", "projection=[Max Intensity]");
	for(k=0; k<iterations; k++) {	
		run("Remove Outliers...", "radius="+radius+" threshold=0 which=Bright");
		//run("Morphological Filters", "operation=Dilation element=Disk radius=3");
	}
	run("Morphological Filters", "operation=Dilation element=Disk radius=10");
	rename("morph");
	selectWindow("MAX_bgstack");
	close();
	selectWindow("morph");
	rename("MAX_bgstack");
	
	for (m=0; m<slices; m++) {
		selectWindow("MAX_bgstack");
		run("Select All");
		run("Copy");
		selectWindow("bgstack");
		setSlice(m+1);
		run("Paste");
	}
	selectWindow("MAX_bgstack");
	close();
	
	imageCalculator("Subtract create stack", "ROFimage","bgstack");
	selectWindow("bgstack");
	close();
	selectWindow("ROFimage");
	close();
	selectWindow("Result of ROFimage");
	rename(imagename);
}	

function ImageFilesOnlyArray (arr) {
	//pass array from getFileList through this e.g. NEWARRAY = ImageFilesOnlyArray(NEWARRAY);
	setOption("ExpandableArrays", true);
	f=0;
	files = newArray;
	for (i = 0; i < arr.length; i++) {
		if(endsWith(arr[i], ".tif") || endsWith(arr[i], ".nd2") || endsWith(arr[i], ".czi") || endsWith(arr[i], ".lsm") ) {   //if it's a tiff image add it to the new array
			files[f] = arr[i];
			f = f+1;
		}
	}
	arr = files;
	arr = Array.sort(arr);
	return arr;
}

function DeleteDir(Dir){
	listDir = getFileList(Dir);
  	//for (j=0; j<listDir.length; j++)
      //print(listDir[j]+": "+File.length(myDir+list[i])+"  "+File. dateLastModified(myDir+list[i]));
 // Delete the files and the directory
	for (j=0; j<listDir.length; j++)
		ok = File.delete(Dir+listDir[j]);
	ok = File.delete(Dir);
	if (File.exists(Dir))
	    print("\\Update10: Unable to delete temporary directory"+ Dir +".");
	else
	    print("\\Update10: Temporary directory "+ Dir +" and files successfully deleted.");
}

function RescaleImage(){
	//Expects FinalRes as an input from user in menu
	input_Title = getTitle();
	input_ID = getImageID();
	//get image information		
	getPixelSize(unit, W, H);
	// Determine rescale value
	Rescale = (1/(FinalRes/W));
	run("Scale...", "x="+Rescale+" y="+Rescale+" interpolation=Bilinear average create");
	rescale_ID = getImageID(); 
	selectImage(input_ID);
	close();
	selectImage(rescale_ID);
	rename(input_Title);
}

function NumberedArray(maxnum) {
	//use to create a numbered array from 1 to maxnum, returns numarr
	//e.g. ChArray = NumberedArray(ChNum);
	numarr = newArray(maxnum);
	for (i=0; i<numarr.length; i++){
		numarr[i] = (i+1);
	}
	return numarr;
}

function closewindow(windowname) {
	if (isOpen(windowname)) { 
      		selectWindow(windowname); 
       		run("Close"); 
  		} 
}



function CloseOpenImages(){
	listwindows = getList("image.titles");
	if (listwindows.length > 0) {
		for (list=0; list<listwindows.length; list++) {
		selectWindow(listwindows[list]);
		close();
		}
	}
}


function CloseOpenWindows(){
	listwindows = getList("window.titles");
	if (listwindows.length > 0) {
		for (list=0; list<listwindows.length; list++) {
			if (listwindows[list] != "Recorder" && listwindows[list] != "Log" && listwindows[list] != "Debug" && listwindows[list] != "B&C") {
				selectWindow(listwindows[list]);
				run("Close");
			}
		}
	}
}


function short_title(imagename){
	nl=lengthOf(imagename);
	nl2=nl-4;
	Sub_Title=substring(imagename,0,nl2);
	return Sub_Title;
}

function CreateSphereAdditive3px(X, Y, Z, Intensity) {
	Radius = 2.5;
	
	run("Colors...", "foreground=white background=black selection=yellow");
	setSlice(Z);
	
	//setColor(Intensity); 
	makeOval(X-Radius+0.5, Y-Radius+0.5, (Radius*2), (Radius*2));
	run("Add...", "value="+Intensity+" slice");
	//print(Intensity);
	
	setSlice(Z+1);
	makeOval(X-Radius+0.5, Y-Radius+0.5, (Radius*2), (Radius*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");
	
	setSlice(Z+2);
	makeOval(X-Radius+1.5, Y-Radius+1.5, ((Radius-1)*2), ((Radius-1)*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");

	setSlice(Z+3);
	makeOval(X-Radius+1.5, Y-Radius+1.5, ((Radius-1)*2), ((Radius-1)*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");

	
	setSlice(Z-1);
	makeOval(X-Radius+0.5, Y-Radius+0.5, (Radius*2), (Radius*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");
	
	setSlice(Z-2);
	makeOval(X-Radius+1.5, Y-Radius+1.5, ((Radius-1)*2), ((Radius-1)*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");

	setSlice(Z-3);
	makeOval(X-Radius+1.5, Y-Radius+1.5, ((Radius-1)*2), ((Radius-1)*2));
	run("Add...", "value="+Intensity+" slice");
	run("Select None");

	
}


function CreateSpot(X, Y, Z, Intensity) {
	run("Colors...", "foreground=white background=black selection=yellow");
	setSlice(Z);
	setColor(Intensity); 
	makeRectangle(X, Y, 1, 1);
	run("Fill", "slice");
}
function append(arr, value) {
	 arr2 = newArray(arr.length+1);
	 for (i=0; i<arr.length; i++)
	    arr2[i] = arr[i];
	 arr2[arr.length] = value;
	 return arr2;
}
