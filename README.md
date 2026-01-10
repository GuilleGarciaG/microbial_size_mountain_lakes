# Resource availability and functional composition underpin size-mediated responses of microbial plankton in mountain lakes
This repository contains codes and data needed to reproduce the analyses and figures of the manuscript:

García-Gómez, G., Sánchez-Hernandez, J., Mas-Gutiérrez, J.A., & Arranz, I. (2026). Resource availability and functional composition underpin size-mediated responses of microbial plankton in mountain lakes.

## Cite the repository as:
García-Gómez, G., Sánchez-Hernandez, J., Mas-Gutiérrez, J.A., & Arranz, I. (2026). Data, code, and statistical analyses of the manuscript: Resource availability and functional composition underpin size-mediated responses of microbial plankton in mountain lakes. (DOI X): [link]

## R scripts
* **picoplankton_data_processing.R** & **nanoplankton_data_processing.R** | These scripts contain all the necessary code to reproduce to process pico- and nanoplankton data, respectively. This workflow obtains raw data from flow cytometry analysis (.fcs files), estimates cell size, and calculates community metrics, including N-M slopes, biovolume, as well as functional dominance variables. The code performs a Principal Component Analysis based on cell autofluorescence, which is used to estimate functional dominance variables. Note that functions to convert light-scatter values from flow cytometry data are demonstrated in the supplementary script **size_calibration.R** (see below).

* **env_data_processing.R** | This script contains code to process and check data of environmental variables, and to perform the Principal Component Analysis (PCA), which is used to quantify resource availability.

* **data_analysis.R** | This script contains code to reproduce statistical models, model diagnostics, and model selection in this study.

* **figures.R** | This script contains code to reproduce figures in the paper, including a brief conceptual description for the hypotheses shown in figure 1.

> ### Supplementary scripts
> **MLE_method.R** | This script contains a custom function to estimate N-M slopes through maximum-likelihood method. It was created by Ignasi Arranz using original functions and materials from package sizeSpectra. See: Edwards, A. M. sizeSpectra (Github, 2019); https://github.com/andrew-edwards/sizeSpectra

> **size_calibration.R** | This script contains code to reproduce the conversion from light-scatter values in flow cytometry data (.fcs files) to cell size (expressed here as "equivalent standard diameter" in µm), using calibration beads of known size.

> **light availability.R** | This script contains code to estimate light availability expressed as light irradiance using open-access satellite data.

## Data
All raw data used in this study are included in folder *raw_data*. Please visit Figshare repository to access the data supporting this study:
[link]

### Microbial plankton
There are two data folders containing all flow cytometry files (.fcs), respectively corresponding to the analysis of picoplankton (*VSSC_26032025_25-473*) and nanoplankton (*FSC_17102024*). These datasets contain .fcs files from the cytometry analysis of each lake sample. 

### Environmental variables
The dataset **environmental_data_lakes_2024.csv** contains measurements of environmental variables for each sample. Columns in this dataset include:

* *ID*: unique ID number for data frame rows
* *massif*: mountain range where sample was collected
* *ID_lake*: unique ID number for each lake
* *ID_sample*: unique ID number for each sample
* *lake*: lake name (in Spanish)
* *habitat*: habitat where sample was collected (littoral: nearshore; pelagic: open water) 
* *replicate*: sample number within each habitat and lake
* *temperature.C*: water temperature, in degrees Celsius
* *dissolved_oxygen.perc*: oxygen saturation, in %
* *total_N_ug.L*: total nitrogen concentration, in µg/l
* *total_P_ug.L*: total phosphorus concentration, in µg/l
* *TOC_mg.L*: total organic carbon, in mg/l
* *nutrilab_label_TN.TP*: ID label from nutrient analyses 
* *nutrilab_label_TOC*: ID label from TOC analyses 

### Lake information
The dataset **lakes_location.csv** contains data on geographic location (latitude and longitude, in decimal degrees) and elevation (in m.a.s.l.) for each lake investigated here. The fine-resolution satellite data of terrain elevation, which was used to estimate light availability, is also available in the same folder ("appRasterSelectAPIService1764233714639-1101222210.tif").

### Size calibration
The folder *calibration_beads* includes two subfolders containing three independent flow cytometry analyses of calibration beads of known diameter, respectively for size-ranges of 0.2-2µm and 1-15µm (with their corresponding folders). Cytometer settings were adjusted to improve measuring precision within each size-range by using violet side-scatter and standard side-scatter light, respectively.

# Notes
All data processing was carried out in the R software version 4.2.2. The *Rcode* folder contains the scripts to reproduce the statistical analyses and figures presented in this manuscript. Please *move the **Rcode** folder inside the global folder (**SI_0126**) before running the code*. Processed data, model outputs, and figures are already stored in their corresponding folders. These already processed and stored data will be used when rerunning statistical analyses and remaking figures.

# R packages
The R packages used for each R script are enlisted in the corresponding R session files (within folder *Rsession*).

# Licence
This repository was provided by the authors under the [to be completed] licence.

# Further information
In case of further questions, please contact: Guillermo García-Gómez, email: guillegar.gz@gmail.com
