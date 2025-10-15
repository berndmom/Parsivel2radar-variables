# Parsivel 2 radar variables

The code is developed to calculate radar variables from measured disdrometer drop size distributions (DSDs) by the Parsivel2. To get from a raw Parsivel2 text file to derived radar variables, three steps are taken. 

## 1. Convert Parsivel raw txt-file to nc-file
The raw Parsivel2 text file is read and converted to a netcdf.

## 2. Filter data and compute drop size distribution parameters
The raw data is taken from Parsivel observations and filtered based on a number of criteria. TBC.
<!-- As Parsivel disdrometers are less sensitive to smaller droplets, size classes with median diameters below 0.4 mm are omitted. The filtering is made to be only applicable to rain and drizzle cases. Therefore, the data are filtered using SYNOP precipitation codes registered by the Parsivel. SYNOP codes are for 51–53 for drizzle, 57-58 for drizzle with rain, and 61–63 for rain, as defined in Table 4680. Additionally, a minimum threshold is set for the detected particles. Furthermore, the data is checked to have at least five consecutive size classes. Size classes falling more than two classes outside this continuity are removed. A closer examination of the remaining data revealed certain features that imply measurement errors, such as insects or side fallers, which deviate in size and fall velocity from typical drops. Raindrops follow a specific velocity-size relationship, as defined by \citet{atlas_doppler_1973}. A limit is set at the 0.5–1.5 percentile of this velocity-size curve for retaining size and velocity classes, with any falling outside being discarded. -->

## 3. Derive radar variables
Radar variables are derived fitting a normalized gamma function to the observed drop size distribution. The found DSD parameters are then used to calculate the corresponding radar variables using the PyTMatrix python packages (https://github.com/jleinonen/pytmatrix).
