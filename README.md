# HDR-US
In this work, we apply HDR techniques to ultrasound imaging, where we combine ultrasound images acquired at different power levels to improve the level of detail visible in the ﬁnal image. Our results strongly suggest that HDR-US imaging can improve the utility of ultrasound in image-based diagnosis and procedure guidance.

### 1. Related Papers
* **High Dynamic Range Ultrasound Imaging**, *A. Degirmenci, D.P. Perrin, and R.D. Howe,* Int J CARS (2018). 
https://link.springer.com/article/10.1007/s11548-018-1729-3 (Free, view-only version: http://rdcu.be/JfOg)

```
@Article{Degirmenci2018,
	author="Degirmenci, Alperen
	and Perrin, Douglas P.
	and Howe, Robert D.",
	title="High dynamic range ultrasound imaging",
	journal="International Journal of Computer Assisted Radiology and Surgery",
	year="2018",
	month="May",
	day="01",
	volume="13",
	number="5",
	pages="721--729",
	issn="1861-6429",
	doi="10.1007/s11548-018-1729-3",
	url="https://doi.org/10.1007/s11548-018-1729-3"
}
```

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
Licensed under the GNU General Public License
Version 3 (GPLv3).
