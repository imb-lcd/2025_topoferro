/*
 *  preprocess_IF.ijm
 *  preporcess the confocol IF images, including scaling, pixel to um and segmentation
 */
 
all_prefix = newArray("YAP-0613-TRULI", "NRF2-0613-KI696-200k");


all_ftype = newArray("c1_sig", "c2_DIC", "c3_dapi");
//all_ftype = newArray("c2_DIC");

for (p = 0; p < lengthOf(all_prefix); p++) {
	for (t = 0; t < lengthOf(all_ftype); t++) {
		
		prefix = all_prefix[p]; //"NRF2-1125-100k-day1";
		ftype = all_ftype[t]; //"c3_dapi";
		print(prefix + " " + ftype);

		File.openSequence("E:/IF/"+prefix+"/"+ftype+"/");
		Stack.setXUnit("um");
		run("Properties...", "channels=1 slices="+nSlices+" frames=1 pixel_width=0.345 pixel_height=0.345 voxel_depth=1.6");
		
		run("Scale...", "x=- y=- z=1.0 width=5000 height=5000 depth=10 interpolation=Bilinear average process create");
		close(ftype);
		selectWindow(ftype+"-1");
		
		File.makeDirectory("E:/IF/" + prefix + "/" + ftype + "_mod/");
		
		run("Image Sequence... ", "select=E:/IF/"+prefix+"/"+ftype+"_mod/ dir=E:/IF/"+prefix+"/"+ftype+"_mod/ format=TIFF use");
		close("*");
	}
}


for (p = 0; p < lengthOf(all_prefix); p++) {

	prefix = all_prefix[p];
	
	File.openSequence("E:/IF/" + prefix + "/c3_dapi_mod/");
	
	File.makeDirectory("E:/IF/" + prefix + "/c3_dapi_density/");
	
	run("Z Project...", "projection=[Max Intensity]");
	
	fname = prefix + "_zMaxc3_ORG.tif";
	
	saveAs("Tiff", "E:/IF/" + prefix + "/c3_dapi_density/" + fname);
	
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+fname+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'Both', 'nTiles':'25', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	
	File.makeDirectory("E:/IF/" + prefix + "/c3_dapi_density/segment/");
	
	stardist_file = prefix + "_zMaxc3_ORG_stardist";
	saveAs("Tiff", "E:/IF/" + prefix + "/c3_dapi_density/segment/" + stardist_file + ".tif");
	
	roiManager("List");
	selectWindow("Overlay Elements of " + stardist_file + ".tif");
	saveAs("Results", "E:/IF/" + prefix + "/c3_dapi_density/" + stardist_file + "_table.txt");
	
	close("ROI Manager");
	close(prefix + "_zMaxc3_ORG_stardist_table.txt");
	close("*");
}