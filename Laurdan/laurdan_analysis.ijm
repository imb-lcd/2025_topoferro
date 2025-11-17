/* 
 * laurdan_analysis.ijm 
 */

/* Modified from the below
 */

/* ImageJ macro for GP images analysis
   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms
   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus
   Nature Protocols 2011
   version July 2011
*/

/* This macro calculates GP images from Laurdan and di-4-ANEPPDHQ ratiometric images in
   bactch mode (whole chosen folder) obtained using a Leica microscope. The generation
   of HSB images of these GP images has been also implemented.
*/

prefix = "laurdan-1019-AY";
lut = "ICA";
Gf = 1;
bc_min = -0.5;	// -1 //	-0.25;	// -1
bc_max = 0.5; // 0.9; //0.7; //	0.3;	// 0.7

setBatchMode(true);

for (w = 9; w <= 9; w++ ) {
//for (w = 2; w <= 2; w++) {
	well = String.pad(w, 2);
	
	print("Well " + well);
	
	// modify input
	dir = "E:/Laurdan/" + prefix + "/Well" + well + File.separator;
	results_dir = "E:/Laurdan/" + prefix + "/Well" + well + File.separator; 
	
	File.makeDirectory(results_dir);
	
	// set up directory
	raw_ordered_images_Dir = dir + "c1_ordered" + File.separator;
	raw_disordered_images_Dir = dir + "c3_disordered" + File.separator;
	
	ordered_images_Dir = dir + "c1_ordered_mod" + File.separator;
	File.makeDirectory(ordered_images_Dir);
	disordered_images_Dir = dir + "c3_disordered_mod" + File.separator;
	File.makeDirectory(disordered_images_Dir);
	
	GP_images_Dir = results_dir + "GP_images" + File.separator;
	File.makeDirectory(GP_images_Dir);
	rawGP_images_Dir = results_dir + "raw_GP_images" + File.separator;
	File.makeDirectory(rawGP_images_Dir);
	mask_images_Dir = results_dir + "Mask_images" + File.separator;
	File.makeDirectory(mask_images_Dir);
	
	HSB_Dir = results_dir + "HSB_images" + File.separator;
	File.makeDirectory(HSB_Dir);

	GP=newArray(256);
	for (i=0; i<256; i++) {
		GP[i]=((i-127)/127);
	}
	
	// Modify the images
	raw_list_ord = getFileList(raw_ordered_images_Dir);
	raw_list_disord = getFileList(raw_disordered_images_Dir);
	
	for (i = 0; i < raw_list_ord.length; i++) {
		open(raw_ordered_images_Dir + raw_list_ord[i]);
		prepareImage("ordered");
		saveAs("tiff", ordered_images_Dir + raw_list_ord[i]);
		close();
		
		open(raw_disordered_images_Dir + raw_list_disord[i]);
		prepareImage("disordered");
		saveAs("tiff", disordered_images_Dir + raw_list_disord[i]);
		close();
	}

	listOrd = getFileList(ordered_images_Dir);
	listDisord = getFileList(disordered_images_Dir);

	for (i = 0; i < listOrd.length; i++) {
//	for (i = 17; i <= 17; i++) {
		chA_name = substring(listOrd[i],0,lengthOf(listOrd[i])-4);
		chB_name = substring(listDisord[i],0,lengthOf(listDisord[i])-4);
	
		// open the ordered and disorded images
		open(ordered_images_Dir + listOrd[i]);
		rename("Image_1a.tif");
		run("Duplicate...","title=Image_1b.tif");
		
		open(disordered_images_Dir + listDisord[i]);
		rename("Image_2a.tif");
		run("Duplicate...","title=Image_2b.tif");
		
		// Perform GP calculation
		// ( Ordered - Disordered ) / (Ordered + Disordered)
		imageCalculator("Substract create 32-bit", "Image_1a.tif", "Image_2a.tif");
		rename("Image_Subs.tif");
		
		imageCalculator("Add create 32-bit", "Image_1b.tif", "Image_2b.tif");
		rename("Image_Add.tif");
		
		imageCalculator("Divide create 32-bit", "Image_Subs.tif", "Image_Add.tif");	
		
		// if divided by 0, set to 0
		changeValues(NaN, NaN, 0);
		
		// Save the raw GP image
		saveAs("tiff", rawGP_images_Dir + chA_name + "_preGP");		
		rename("Image_preGP.tif");
			
		// Set the range to the GP value range
		// then automatically scale for the display range of 16-bit images
		setMinAndMax(-1.0000, 1.0000);
		call("ij.ImagePlus.setDefault16bitRange", 0);
		
		// Generate 1 bit mask
		selectImage("Image_Add.tif");
		run("Duplicate...","title=Image_1bit.tif");
//		setThreshold(t*2, 510); // original code, t = 50. 
		setAutoThreshold("Default dark");
		run("Convert to Mask");
		run("Subtract...", "value=254");
		saveAs("tiff", mask_images_Dir + chA_name + "_1bitmask");
		rename("Image_1bit.tif");

		// Generate GP image that removes the 1 bit mask
		imageCalculator("Multiply create 32-bit", "Image_1bit.tif", "Image_preGP.tif");

//		setMinAndMax(-0.2, 0.60000);

		run(lut);
		saveAs("tiff", GP_images_Dir + chA_name + "_GP");
		
		//
		// Generate HSV image
		//
		// Make the Intensity channel (Ordered+Disordered) as "brigthness"
		selectImage("Image_Add.tif");
		run("Enhance Contrast", "saturated=0.5 normalize");
		rename("Brightness");

		// Make the preGP image as "Hue" and adjust B/C
		selectImage("Image_preGP.tif");
		run(lut);
		rename("Hue");
		
		run("Brightness/Contrast...");
		setMinAndMax(bc_min, bc_max);
		
		// Split the image by RGB, and multiply the intensity channel
		run("RGB Color");
		run("Split Channels");

		imageCalculator("Multiply create 32-bit", "Brightness", "Hue (red)");
		rename("bR");
		run("8-bit");

		imageCalculator("Multiply create 32-bit", "Brightness", "Hue (green)");
		rename("bG");
		run("8-bit");

		imageCalculator("Multiply create 32-bit", "Brightness", "Hue (blue)");
		rename("bB");
		run("8-bit");

		// merge the RGB channels
		run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
		saveAs("tiff", HSB_Dir + chA_name + "_HSB");		

		closeAllImages();
	}
}

///////////////// FUNCTIONS ////////////////////

function closeAllImages() {				// This function closes all images
	while (nImages>0) {
		selectImage(nImages);
		close();
	}
}

function newFolder() {					// This function creates a folder, removing any existing file in a folder with the same name
	File.makeDirectory(Folder);
	listFolder = getFileList(Folder);
	for (i = 0; i < listFolder.length; i++) {
		File.delete(Folder+listFolder[i]);
	}
}

function prepareImage (order_type) {				// This funcion prepares each image for the analysis
	s=getTitle;
//	run("8-bit");
	run("Grays");
	run("32-bit");
	run("Median...", "radius=1");

	return s;
}
