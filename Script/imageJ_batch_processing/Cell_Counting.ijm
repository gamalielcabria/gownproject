// Defining Functions:
// Returns true if the window with the given title is open; otherwise false.
function isWindowOpen(title) {
    // Attempt to select the window (this will fail if the window does not exist)
    // exec() returns an empty string if successful, or an error message if not.
    err = exec("selectWindow(\"" + title + "\");");
    if (err == "") {
        // We successfully selected it, so the window is open.
        return true;
    } else {
        // The window doesn't exist or can't be selected.
        return false;
    }
}

// Define the names for each box (for saving filenames)
boxes = newArray("Box1", "Box2", "Box3", "Box4", "Box5", "Box6");

// Define the coordinates for each rectangle:
//   X, Y, Width, Height
// (ImageJ macros do not have a direct 2D array type; we store each dimension in a separate 1D array.)
roisX = newArray(   33, 39, 924,  933,  1821, 1833);
roisY = newArray(   93, 978, 75,  966,    72,  960);
roisW = newArray(   894, 897, 912,  900,   876,  867);
roisH = newArray(   879, 891, 891,  894,   891,  891);

// Define input and output directories when running this batch run within ImageJ 
inputDir = getDirectory("Choose Input Directory");

// Get the list of files in the input directory
fileList = getFileList(inputDir);

// Loop through each file in the input directory
for (j = 0; j < fileList.length; j++) {
    // Get the current file name
    fileName = fileList[j];

    if (endsWith(fileName, "_edited.tif")) {
        // Get the Prefix of the filename
        imageName = replace(fileName, "_edited.tif", "");

        // Open the current image
        open(inputDir + fileName);
        originalImageTitle = getTitle();

        // Set the image scale and grayscale
        run("8-bit");
        run("Set Scale...", "distance=35642 known=10000 unit=micron");

        // Loop through each box
        for (i = 0; i < boxes.length; i++) {
            // Re-select the original image each time
            selectImage(originalImageTitle);
			roiManager("reset");            
            run("Duplicate...", "title=" + imageName + "_Box" + (i+1) );
            
            // Make the rectangle for the i-th ROI
            makeRectangle(roisX[i], roisY[i], roisW[i], roisH[i]);

            // Analyze particles
            run("Analyze Particles...", "size=1-Infinity circularity=0.80-1.00 display exclude clear summarize overlay add");

            // Save the overlay for the current box
            // overlay = Overlay.copy();
            // roiManager("Add");
            
            // Save the Results table
            particleCount=roiManager("count");
            // err = eval("selectWindow(\"Results\");");
            print("Results Check: " +  " error-end \nParticle count: " + particleCount);
            if (particleCount > 0 ){
                selectWindow("Results");
                saveAs("Results", inputDir + imageName + "_Box" + (i+1) + "_results.csv");
                run("Close");
            }

            // Apply overlay options and flatten the image
            run("Overlay Options...", "stroke=red width=1 fill=red set apply");
            run("Flatten");

            // Save the overlay image
            saveAs("Tiff", inputDir + imageName + "_Box" + (i+1) + "_overlay.tif");
        }

        // Save the Summary table per image
        err = eval("macro","selectWindow(\"Summary\");");
        print("Summary Check: " + err + " error-end");
        if (err == "0" ){
            selectWindow("Summary");
            saveAs("Results", inputDir + imageName + "_count.csv");
            run("Close");
        }

        // Close the current image
        run("Close All");
    }
}