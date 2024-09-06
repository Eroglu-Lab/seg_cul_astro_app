# seg_cul_astro_app

## How to run seg_astro.mlapp
1. Install the latest version of MATLAB
1. Go to the folder ‘seg_cul_astrocyte_app’ and double click ‘seg_astro.mlapp’
1. A new window will pop-out, the first row of the left panel is for input the image data (Input the two channels separately), the second row and the third row are for the detection. 
[insert image 1]
1. After loading the data, a folder will be automatically created along side the image data, where the resulted will be stored.
[insert image 2]
[insert image 3]
1. Then press the buttons for ‘detect mitochondria’ and ‘detect astrocyte’. Detecting astrocyte can take around 10 mins, depending on the configuration of your PC.	
[insert image 4]
[insert image 5]
1. You can also check the results in the folder
[insert image 6]
1. You can also skip the two detection steps by input the mask of the astrocyte ‘astro_detection_ssp.tif’  as well as the mask of the mitochondria ‘mito_detection.tif’and the image of the astrocyte (if the detection has already been done and you want to redo the branch detection).
[insert image 7]
1. Once the branch detection is done (default parameter is 1.5), you can change the width threshold. First input a number and then press Enter.
[insert image 8]
1. The next step is to refine the branches (some main branches might be isolated; this step is to recover some isolated branches and link them back to the core branches)
1. Newly added functions: detecting the second layer branches. After detecting the root branches, the secondary branches are again searched.
[insert image 9] 
the red region belongs to the secondary branches, while the blue region belongs to the tertiary branches. The green region is the detected wide part in the previous step.
1. The final refine branches step is to remove any distant/ remote branches
	[insert image 10]
1. The quantification step analyzes the number/ size/ circularity of mitochondria at each level of branches, which are primary, secondary, fine branches and terminal branches. After press ‘run’ button, eight files will be generated at the same folder where the input is located.  
[insert image 11]
[insert image 12]
1. Four of the files are the summary of the mitochondria in different types of branches, the name starts with ‘mito’. This .txt file can also be opened in Excel. The number/ mean size/ mean circularity is displayed in the first three rows. And starting from fifth row, the size and circularity of each mitochondrion is calculated. Circularity: 4πA/c^2 , where A is the size of the area and c is the perimeter of the region. The other four .txt files are the summary of the branches, containing the area and length of each type of branch.
[insert image 13]
