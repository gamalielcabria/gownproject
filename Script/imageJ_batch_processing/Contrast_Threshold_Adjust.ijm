// Running the Macro as a headless Script
// Example: Fiji.app/ImageJ-linux64 --headless -macro process_images.ijm "/input/path,/output/path"
// args = getArgument();
// splitArgs = split(args, ",");
// inputDir = splitArgs[0];
// outputDir = splitArgs[1];

// Define input and output directories when running this batch run within ImageJ 
inputDir = getDirectory("Choose Input Directory");
outputDir = getDirectory("Select Output Directory");

// Get the list of files in the input directory
fileList = getFileList(inputDir);

// Loop through each file in the input directory
for (i = 0; i < fileList.length; i++) {
    // Get the current file name
    fileName = fileList[i];
    
    // Check if the file is a .jpg or .jpeg or .png
    if (endsWith(fileName, ".jpg") || endsWith(fileName, ".jpeg" ) || endsWith(fileName, ".png" )) {
        // Opening the File
        open(inputDir + fileName);

        // Main Macros for ImageJ
        run("Duplicate...", "title=" + fileName + "_edited");
        selectImage(fileName);
        close;
        selectImage(fileName + "_edited");
        run("8-bit");
        run("Median...", "radius=2");
        run("Enhance Contrast...", "saturated=0.30");
        run("Subtract Background...", "rolling=50 sliding");
        run("Enhance Contrast...", "saturated=0.50");
        if (fileName == "D2-SMPL-DAPI-BR1-10X-1.png") {
            //run("Enhance Contrast...", "saturated=2.50");
            setAutoThreshold("Default dark");
            //run("Threshold...");
            setThreshold(54, 255);
        } else {
            setAutoThreshold("Default dark");
            //run("Threshold...");
            setThreshold(19, 255);
        }
        setOption("BlackBackground", true);
        run("Convert to Mask");
        //run("Close");
        run("Dilate");
        //run("Watershed");
        
        // Saving edited File
        if (endsWith(fileName, ".jpg")) {
            newName = replace(fileName, ".jpg", "_edited.tif");
        } else if (endsWith(fileName, ".jpeg")) {
            newName = replace(fileName, ".jpeg", "_edited.tif");
        } else if (endsWith(fileName, ".png" )) {
            newName = replace(fileName, ".png", "_edited.tif");
        }
        
        // Save the processed image as a new .tif file
        saveAs("Tiff", outputDir + newName);
        
        // Close the images
        close();
    }
}
