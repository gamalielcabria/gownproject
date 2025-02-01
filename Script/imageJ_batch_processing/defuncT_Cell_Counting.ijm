
run("8-bit")
run("Set Scale...", "distance=35642 known=10000 unit=micron");
run("Overlay Options...", "set");
run("ROI Manager...");
roiManager("Reset");

// Upper Left Box - Box1
makeRectangle(33, 93, 894, 879);
run("Analyze Particles...", "size=1-Infinity circularity=0.80-1.00 display exclude clear summarize overlay add");
firstOverlay = Overlay.copy();
roiManager("Add");

selectWindow("Results");
saveAs("Results", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box1_results.csv");
run("Close");

selectWindow("Summary");
saveAs("Results", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box1_count.csv");
run("Close");

run("Overlay Options...", "stroke=red width=1 fill=red set apply");
run("Flatten");
saveAs("Tiff", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box1_overlay.tif");


// Lower Bottom Box - Box2
makeRectangle(39, 978, 897, 891);
run("Analyze Particles...", "size=1-Infinity circularity=0.80-1.00 display exclude clear summarize overlay add");
secondOverlay = Overlay.copy();
roiManager("Add");

selectWindow("Results");
saveAs("Results", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box2_results.csv");
run("Close");

selectWindow("Summary");
saveAs("Results", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box2_count.csv");
run("Close");

run("Overlay Options...", "stroke=red width=1 fill=red set apply");
run("Flatten");
saveAs("Tiff", "/home/glbcabria/Workbench/P3/Results/CellcountPhotos/Training/D2-Edited_Training/D2-SMPL-DAPI-GR1-10X-3_Box2_overlay.tif");

// Upper Center Box - Box3
//makeRectangle(924, 75, 912, 891);
// Lower Center Box - Box4
//makeRectangle(933, 966, 900, 894);
// Upper Right Box - Box5
//makeRectangle(1821, 72, 876, 891);
// Lower Right Box - Box6
//makeRectangle(1833, 960, 867, 891);

