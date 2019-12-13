# scripts for saving & plotting region-of-interest (ROI) timecourses for fmri data

These matlab scripts provide an example of how to save out & then plot ROI timecourses. It includes sample data from 3 subjects on the MID task. 

...turns out that the sample fmri data exceeds the file size limit for git! Please email Kelly at hennigan@stanford.edu if you'd like access to the sample fmri data. 

## Getting started



### Software requirements 

* [Matlab](https://www.mathworks.com/products/matlab.html)
* [matlab package, VISTASOFT](https://github.com/vistalab/vistasoft) (their niftiRead() function for loading nifti data is required)


## Pipeline 

In matlab, run:
```
saveOutRoiTimeCourses_script
```
to save out timecourses, averaged over trials for each condition.

then run:
```
plotRoiTimeCourses_script
```
to plot timecourses (averaged across trials). 





