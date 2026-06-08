// ── Constants ──────────────────────────────────────────────────────────────
WORK_DIR    = "T:/Jen-Hao/Laurdan/";
SCRIPT_DIR  = WORK_DIR + "scripts_ccwu19/";
RIDGE_DIR   = SCRIPT_DIR + "ridges/";
INPUT_CSV   = SCRIPT_DIR + "fname_slices_thres.csv";

RIDGE_PARAMS = "line_width=29 high_contrast=300 low_contrast=100 " +
               "extend_line displayresults add_to_manager make_binary " +
               "method_for_overlap_resolution=NONE sigma=8.87 " +
               "lower_threshold=0 upper_threshold=0.17 " +
               "minimum_line_length=100 maximum=0";

// ── 0. Cleanup ──────────────────────────────────────────────────────────────
close("Results");
close("ROI Manager");
close("*.csv");
close("*");

// ── 1. Load input table ─────────────────────────────────────────────────────
open(INPUT_CSV);
nRows = Table.size;
//setBatchMode("hide");

// ── 2. Main loop ────────────────────────────────────────────────────────────
for (i = 0; i < nRows; i++) {
    print("Processing [" + (i+1) + "/" + nRows + "]");

    // -- Read row metadata --
    selectWindow("fname_slices_thres.csv");
    thres  = Table.get("thres",    i);
    sliceN = Table.get("slice",    i);   // numeric slice index
    slice  = IJ.pad(sliceN, 2);          // zero-padded string for filenames
    inPath = Table.getString("location", i);

    // -- Open image and select slice --
    open(inPath);
    setSlice(sliceN);

    // -- Threshold → binary mask --
    run("Duplicate...", "use");
    setAutoThreshold("Default dark");
    setThreshold(thres, 1e30);
    setOption("BlackBackground", true);
    run("Convert to Mask");

    // -- Ridge detection --
    run("Ridge Detection", RIDGE_PARAMS);

    // -- Discard intermediate windows --
    close("Results");
    close("Junctions");
    close("Summary");
    close("*");

    // -- Save ROI zip --
    oname    = File.getName(inPath);
    idx_All  = indexOf(oname, "All");
    roi_file = substring(oname, 0, idx_All) + slice + "_ridge.zip";
    roiManager("Save", RIDGE_DIR + roi_file);
    close("ROI Manager");
}

print("Done.");