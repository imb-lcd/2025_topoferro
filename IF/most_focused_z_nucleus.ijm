/*
 * identify_focused_z_by_nucleus
 * identify the most focused z by the nucleus signals
 */

ROOTDIR = "E:/IF/";

focus_type = "sig"; // "c1_sig"

if (focus_type == "dapi") {
	counts = 5;
} else {
	counts = 4;
}

all_prefix = newArray("YAP-0613-TRULI", "NRF2-0613-KI696-200k");


//start_z = newArray(4, 3, 4, 4, 4, 4, 4);
//end_z   = newArray(7, 6, 7, 7, 7, 7, 7);
start_z = newArray(1, 1);
all_end_frame = newArray(15, 16);;


for(p = 0; p < lengthOf(all_prefix); p++) {
	prefix = all_prefix[p]; //"NRF2-1125-100k-day3";
	
	// load IF path
	
	if (focus_type == "sig") {
		IFPATH = "c1_" + focus_type + "_mod";
	} else {
		IFPATH = "c3_" + focus_type + "_mod";
	}
	
	if_file =  ROOTDIR + prefix + "/" + IFPATH + "/";
	File.openSequence(if_file);
	run("16-bit");
	
	// load nucleus locations
	DENPATH = "c3_dapi_density/";
	File.makeDirectory(ROOTDIR + prefix + "/" + DENPATH + "/");
	file = ROOTDIR + prefix + "/" + DENPATH + "/" + prefix + "_zMaxc3_ORG_stardist_table.txt";
	fileContent = File.openAsString(file); 
	lines = split(fileContent, "\n"); 
	
	lineLength = lengthOf(lines);
	
	outfile = ROOTDIR + prefix + "/" + DENPATH + "/" + prefix+ "_zMaxc3_ORG_stardist_table_focus_" + focus_type + ".txt";

	out = File.open(outfile);
	
	header = "index\tName\tType\tGroup\tX\tY\tWidth\tHeight\tPoints\tColor\tFill\tLwidth\tPos\tC\tZ\tT\tFocus";
	print(out, header);

	
	// for every single cell
	for (i = 1; i < lineLength; i++) {
		line = split(lines[i], "\t");
		
		x = parseFloat(line[5])-10;
		y = parseFloat(line[6])-10;
		width = parseFloat(line[7])+10;
		height = parseFloat(line[8])+10;
		
		x = Math.max(x, 0);
		y = Math.max(y, 0);
	
		// select IF window and region
		selectWindow(IFPATH);
		makeRectangle(x, y, width, height);
		
		print(x, y, width, height);
		
		// calculate the most focused slice
		sz = 9;
		run( "Focus LP", "criterion=Standard size=" + sz);
		//Macro code for accessing the result:
		returnStr = split( call( "Focus_LP.macroReturn" ), " " );
		focusedSlice = parseInt( returnStr[0] );
	
		outline = newArray(lengthOf(line) - 1);
		for(k = 0; k < lengthOf(line)-1; k++) {
			outline[k] = line[k+1];
		}
	
		outline = String.join(outline, "\t")+"\t"+(focusedSlice+start_z[p]-1);
	
		print(out, outline);
	}
	File.close(out);
	close("*");
}
