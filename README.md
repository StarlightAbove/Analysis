# Analysis

A repository containing all code and data used for the project started over the summer at the Przybyl Lab, including all data, plots and code.

The final code should be in the "Comparative Code" folder, the rest of the code is experimental. fpCheck is for LM data, and gRangeAccModel is for LMS models. 

visualize.R is the actual visualizing code.

# Indexes for all statistics
### "LMS" - LMS cases from TCGA-provided samples.
| Sample Number (Internal) | TCGA Identifier | File Path                                    | Most/Least Accurate? | Comments |
|--------------------------|-----------------|----------------------------------------------|----------------------|----------|
| Sample 1                 | TCGA-3B-A9HQ    | ./Statistics/LMS/AllSamples/TCGA-3B-A9HQ.png |                      |          |
| Sample 2                 | TCGA-3B-A9HX    | ./Statistics/LMS/AllSamples/TCGA-3B-A9HX.png |                      |          |
| Sample 3                 | TCGA-DX-A3UC    | ./Statistics/LMS/AllSamples/TCGA-DX-A3UC.png |                      |          |
| Sample 4                 | TCGA-DX-A48O    | ./Statistics/LMS/AllSamples/TCGA-DX-A48O.png |                      |          |
| Sample 5                 | TCGA-DX-A48P    | ./Statistics/LMS/AllSamples/TCGA-DX-A48P.png |                      |          |
| Sample 6                 | TCGA-DX-A6BA    | ./Statistics/LMS/AllSamples/TCGA-DX-A6BA.png | ✖ Least Accurate     | Labelled as lAcc.png. |
| Sample 7                 | TCGA-DX-A6Z2    | ./Statistics/LMS/AllSamples/TCGA-DX-A6Z2.png |                      |          |
| Sample 8                 | TCGA-DX-A7EL    | ./Statistics/LMS/AllSamples/TCGA-DX-A7EL.png |                      |          |
| Sample 9                 | TCGA-DX-A7EN    | ./Statistics/LMS/AllSamples/TCGA-DX-A7EN.png | ✔ Most Accurate      | Labelled as mAcc.png  |
| Sample 10                | TCGA-K1-A42W    | ./Statistics/LMS/AllSamples/TCGA-K1-A42W.png |                      |          |

### "LabLMS" - LMS cases from laboratory-provided samples.
| Sample Number (Internal) | STT Identifier | Filepath to Manhattan Plot              | Most/Least Accurate | Comments |
|--------------------------|----------------|-----------------------------------------|---------------------|----------|
| Sample 1                 | 9202           | ./Statistics/LabLMS/AllSamples/9202.png |                     |          |
| Sample 2                 | 9203           | ./Statistics/LabLMS/AllSamples/9203.png | ✔ Most Accurate     |          |
| Sample 3                 | 9327           | ./Statistics/LabLMS/AllSamples/9327.png |                     |          |
| Sample 4                 | 9328           | ./Statistics/LabLMS/AllSamples/9328.png |                     |          |
| Sample 5                 | 9337           | ./Statistics/LabLMS/AllSamples/9337.png |                     |          |
| Sample 6                 | 9338           | ./Statistics/LabLMS/AllSamples/9338.png |                     |          |
| Sample 7                 | 9350           | ./Statistics/LabLMS/AllSamples/9350.png |                     |          |
| Sample 8                 | 9353           | ./Statistics/LabLMS/AllSamples/9353.png | ✖ Least Accurate    |          |
| Sample 9                 | 9354           | ./Statistics/LabLMS/AllSamples/9354.png |                     |          |
| Sample 10                | 9355           | ./Statistics/LabLMS/AllSamples/9355.png |                     |          |
| Sample 11                | 9356           | ./Statistics/LabLMS/AllSamples/9356.png |                     |          |
| Sample 12                | 9357           | ./Statistics/LabLMS/AllSamples/9357.png |                     |          |
| Sample 13                | 9358           | ./Statistics/LabLMS/AllSamples/9358.png |                     |          |

### "LM" - LM cases from laboratory-provided samples.
| Sample Number (Internal) | STT Identifier      | Filepath to Manhattan Plot                         | Most/Least Accurate | Comments |
|--------------------------|---------------------|----------------------------------------------------|---------------------|----------|
| Sample 1                 | 203219650162_R08C01 | ./Statistics/LM/AllSamples/203219650162_R08C01.png |                     |          |
| Sample 2                 | 203219650188_R08C01 | ./Statistics/LM/AllSamples/203219650188_R08C01.png | ✖ Least Accurate    |          |
| Sample 3                 | 203219750019_R05C01 | ./Statistics/LM/AllSamples/203219750019_R05C01.png |                     |          |
| Sample 4                 | 203219750053_R05C01 | ./Statistics/LM/AllSamples/203219750053_R05C01.png |                     |          |
| Sample 5                 | 203219760009_R08C01 | ./Statistics/LM/AllSamples/203219760009_R08C01.png |                     |          |
| Sample 6                 | 203220070150_R02C01 | ./Statistics/LM/AllSamples/203220070150_R02C01.png |                     |          |
| Sample 7                 | 203225140266_R03C01 | ./Statistics/LM/AllSamples/203225140266_R03C01.png |                     |          |
| Sample 8                 | 203257030015_R06C01 | ./Statistics/LM/AllSamples/203257030015_R06C01.png |                     |          |
| Sample 9                 | 203257030060_R08C01 | ./Statistics/LM/AllSamples/203257030060_R08C01.png | ✔ Most Accurate     |          |

Most accurate Manhattan plots are usually labelled as mAcc.png, and least accurate as lAcc.png.
