# 2D_peakshape
GUI for fitting of sum of multivariate normal distributions to series of 2D peaks


![gif](https://github.com/nowacowski/2D_peakshape/blob/master/usage.gif)
<br>

- Fit2D
  - example of usage, by creating series of two overlapping peaks
  - passing to prepare_input.m to make sure input for Fit2D_gui.m is correct
  - fitting multivariate gaussians using created GUI Fit2D_gui.m
  - creating model from fitting parameters using Gauss2D.m
  
- prepare_input
  - taking 3D array with amplitude data, and vectors with values of x, y and z axes
  - making sure all data is provided and correct and transform if needed for proper input for Fit2D_gui.m
  
- Gauss2D
  - function creating sum of multivariate normal distributions
  - requires input with axes grid and parameters for multivariate Gaussian
  
- caruana
  - function to analytically fit Gaussian to data and output obtained parameters
  - based on Caruana's algorithm <br>
  *R. Caruana, R. Searle, T. Heller, and S. Shupack, “Fast algorithm for the resolution of spectra,” Anal. Chem., vol.58, no. 6, pp. 1162–1167, May 1986*
  
- cov_error
  - function to calculate error obtained from the fitting using Jacobian obtained from lsqcurvefit function
  
- Fit2D_gui
  - GUI function for interactive fitting and visualizing of fitted data
  - usage showed on gif (usage.gif) and described within script
  
- InParamGui
  - GUI function for choosing initial parameters for fitting
  - used within Fit2D_gui as a separate window
  - empty table or filled with previously input initial parameters or latest result from fitting of first slice from data series
