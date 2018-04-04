# HDR-US
High Dynamic Range Ultrasound Imaging

### 1. Related Papers
* **High Dynamic Range Ultrasound Imaging**, *A. Degirmenci, D.P. Perrin, and R.D. Howe,* Int J CARS (2018). 
https://link.springer.com/article/10.1007/s11548-018-1729-3 (Free, view-only version: http://rdcu.be/JfOg)

### 2. Installation
Download the repo to your machine:

	git clone https://github.com/adegirmenci/HDR-US.git
	
From the root project directory, run the install script:

	installHDRUS
  
#### 2.1 MATLAB Toolbox Dependencies
* Image Processing Toolbox
* Parallel Computing Toolbox (optional)

This code was tested using MATLAB 2016b and 2017b.

#### 2.2 Dependencies

We use F. Banterle's HDR Toolbox (https://github.com/banterle/HDR_Toolbox). The DebevecCRF function is modified to return two extra variables, logE and stack samples.

### 3 Usage
Run the MATLAB script runHDRUS.m:

	runHDRUS

Select the dataset directory from the dialog.

To save results, set the flag

	saveResults = true;

### 5 License
HDR-US was developed at the Harvard Biorobotics Lab.
The  licensed under the GNU General Public License
Version 3 (GPLv3).
