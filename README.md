# Magnetic-Dipole-Spatial-Coordinate-Inversion-Algorithm
## Statement 
This code demonstrates an inversion algorithm for solving the position and magnetic moment parameters of a dipole through coordinate transformation. It serves the paper titled "Magnetic Dipole Spatial Coordinate Transformation and Inversion Method" with Lu Yin as the first author, and is shared here for communication and learning purposes.
## Software Download and Installation 
1. The software used in this work is MATLAB (R2021b, 64-bit), running on Windows 10. The MATLAB license is owned by the University of Science and Technology Beijing (USTB), the author’s affiliated institution.
2. Note for Users Without MATLAB Access:
If you do not have a MATLAB license, please install it through official channels ，or refer to the pseudocode provided (Pseudocode.txt) in the supplementary files to implement the algorithms in other programming languages.
## Code Usage Instructions
### Parameter Configuration Background
In the example, the cuboid magnetic source has center coordinates at (10, -7, -1) with dimensions of 6.5 × 4 × 20. Its xy cross-sectional view is shown in the figure：
![image](https://github.com/user-attachments/assets/7b8d308c-cb1e-4123-833f-686a0d9d567e)
### Code Usage Instructions
1. For users who have MATLAB installed, simply run the main script "text_example.m". 
P.S. 1）Ensure that the main file and all supporting function files are in the same folder during execution.
     2）Before running the file, please modify the data reading path for Bz (Line 16 in the code).
2. After running the script, two figures will be displayed directly:
       Figure 1 shows the loaded Bz data.
       Figure 2 presents the inversion results.
   ![image](https://github.com/user-attachments/assets/eb982d77-e65d-4b86-80fa-5262462cbe98)

3. The Fourier-based gradient calculation is implemented in the function file "text_Fourier.m", while the inversion computation is performed in the function file "text_inverse.m".
## If you encounter any issues
Please contact me ：yinlu_1999@163.com
## License
The author's affiliated institution, the University of Science and Technology Beijing (USTB), holds a licensed MATLAB campus-wide agreement. For details, please refer to: https://soft.ustb.edu.cn/product.html?id=320
