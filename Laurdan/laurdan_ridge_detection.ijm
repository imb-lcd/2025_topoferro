/* 
 * Quantify lipid-orderedness through ridge detection
 * laurdan_ridge_detection.ijm
 */
 
prefix = "laurdan-0706-control";

setColor("yellow");
setLineWidth(3);


w_start = 3;
w_end = 3;

welllist = newArray(7, 13);

dir = "E:/Laurdan/";
outfile = dir + prefix + "_ridge.txt";

out = File.open(outfile);

for (w = w_start; w <= w_end; w++) {
	well = String.pad(w, 2);
	print(well);
	
	im_path = dir + prefix + "_crop_cell/Well" + well + "/";

	im_list = getFileList(im_path);
	
//	for (i = 0; i < im_list.length; i++) {
	for (i = 10; i < 20; i++) {
		fname = im_list[i];
		print(fname);
		
		open(im_path + "/" + fname);
		
		// run ridge detection
		run("Ridge Detection", "line_width=29 high_contrast=300 low_contrast=100 extend_line displayresults add_to_manager method_for_overlap_resolution=NONE sigma=8.87 lower_threshold=0 upper_threshold=0.17 minimum_line_length=75 maximum=0");

		// run ridge detection
//		run("Ridge Detection", "line_width=20 high_contrast=150 low_contrast=50 extend_line displayresults add_to_manager method_for_overlap_resolution=NONE sigma=6.27 lower_threshold=0 upper_threshold=0.17 minimum_line_length=75 maximum=0");
		
		close("Junctions");
		close("Results");
		close("ROI Manager");
	
		Table.rename("Summary", "Results");
	
		len = 0;
		for (row = 0; row < nResults; row++) {
			len += getResult("Length", row);
		}
		
		outline = "";
		idx = indexOf(fname, "cell")+4;
		cell = substring(fname, idx, idx+2);

		outline = fname + "\t" + well + "\t" + cell + "\t" + len;
		print(out, outline);
	}
}
File.close(out);
