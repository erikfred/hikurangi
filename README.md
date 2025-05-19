HIKURANGI
=========

These codes ingest Hikurangi seafloor pressure data into Matlab and process the data with decimation routines, tidal filters, and options for drift- and seasonal-signal removal. There is not a lot of organization here, as a clear project did not emerge prior to switching gears to work of A-0-A drift-corrected data.

Using on your local machine
---------------------------

Throughout, the paths are somewhat generalized (e.g., '../figures/blahblahblah.png), but in practice you'll likely get directory errors unless you manage to duplicate my unkown-to-you directory structure. So, running these codes will require you to update the paths to suit your system and organizational preferences.

THE EXCEPTION to this are the data import codes, which I have tried to make very user friendly with path variables defined at the beginning of the codes.

Code descriptions
----------

* /apg_read

*(Originally authored by Neville Palmer) Directory with Python scripts for converting GNS APG binary files to text files*

* /old

*Directory with miscelaneous codes that were used for a specific purpose that has passed, or that have been superceded by other codes*

* /tools <--- add to the path to ensure main codes work!

*Useful bits of code (many authored by others) that get called by the main codes*

* CPIES_maps.m

*Use ECCO2 ocean model output and CPIES locations to explore which might be the most useful in the context of the A-0-A/secular strain problem*

* depth_dist_RMS.m

*For the various experiment years, plot depth difference vs distance vs RMS of differenced pressures for all possible station pairs. Goal is to assess how effective depth-matched differencing is and how sensitive it is to changing depths and increasing separation between stations*

* detect_ramps.m (INCOMPLETE)

*Attempting to apply ramp detector from Fredrickson et al. (2023) to Hikurangi data*

* difference_all.m

*For the various experiment years, difference all possible station combinations and make plots. Creates fodder for visually searching for interesting transient signals*

* get_ts.csh

*Shell script for pulling GNSS data from GNS servers*

* H7_offsets.m

*Inspects a few HOBITSS VII (2020-2021) pressure records at finer detail to try to understand some spurious offsets seen in the data*

* import_HX.m (or _GX.m)

*Reads the original data files (usually some kind of text file) for each instrument in the experiment year. Observations are trimmed to take out the deployment and recovery intervals. When necessary, pressure data are manipulated to correct anomalies such as offsets. The cleaned up data are saved in several formats and plotted.*

* method_comp_plots.m

*For HOBITSS VIII experiment (2021-2022), compares depth matched differences vs network average differences vs CEOF corrections and makes plots*

* seasonal_correction.m

*For any input pressure data, applies four flavors of drift/seasonal signal correction: 1) 3rd-order polynomial + exponential, 2) sinusoid + exponential, 3) linear + exponential, 4) CEOF + exponential*

* simple_maps.m

*For the various experiment years, generates map-view plots of station distributions*