#! /bin/bash

# No useful data recorsed as incorrect APG cable was installed.
# python ~/scripts/seafloor-geod/apg_read.py \
#     --infile        ./TAN2301_GNS22-PF_140338.apg \
#     --snapg         140338 \
#     --version       CSAC2013 \
#     --clkstart      2022-10-15_00:20:00 \
#     --endwndw       2022-10-16_00:00:00.0 \
#     --bininterval   1d \
#     --gpssynctime   2023-006:12:15:00 \
#     --synctickcount 0x01adfe6e57 \
#     --fmttime       d \
#     | tee ./out/GNS22-PF.log
    # --tempsmth      5001 \
    # --decimate      1 \
    # --beginwndw     2022-10-16_00:00:00.0 \
    # --plot          rs \
    # --outfile       ./out/GNS22-PF.csv \

python ~/scripts/seafloor-geod/apg_read.py \
    --infile        ./TAN2301_GNS22-PI_140339.apg \
    --snapg         140339 \
    --version       Seascan2018 \
    --tempsmth      5001 \
    --decimate      1 \
    --clkstart      2022-10-13_23:29:00 \
    --bininterval   1d \
    --gpssynctime   2023:006:06:49:00 \
    --fmttime       d \
    --outfile       ./out/GNS22-PI.csv \
    | tee ./out/GNS22-PI.log
    # --beginwndw     2022-11-01_00:00:00.0 \
    # --endwndw       2022-11-02_00:00:00.0 \
    # --plot          rs \

python ~/scripts/seafloor-geod/apg_read.py \
    --infile        ./TAN2301_GNS22-PJ_140345.apg \
    --snapg         140345 \
    --version       Seascan2018 \
    --tempsmth      5001 \
    --decimate      1 \
    --clkstart      2022-10-12_23:28:30 \
    --bininterval   1d \
    --gpssynctime   2023:006:03:01:00 \
    --fmttime       d \
    --outfile       ./out/GNS22-PJ.csv \
    | tee ./out/GNS22-PJ.log
    # --beginwndw     2022-11-01_00:00:00.0 \
    # --endwndw       2022-11-02_00:00:00.0 \
    # --plot          rs \

python ~/scripts/seafloor-geod/apg_read.py \
    --infile        ./TAN2301_GNS22-PK_140344.apg \
    --snapg         140344 \
    --version       Seascan2018 \
    --tempsmth      5001 \
    --decimate      1 \
    --clkstart      2022-10-11_23:28:00 \
    --bininterval   1d \
    --gpssynctime   2023:007:01:48:30 \
    --fmttime       d \
    --outfile       ./out/GNS22-PK.csv \
    | tee ./out/GNS22-PK.log
    # --beginwndw     2022-11-01_00:00:00.0 \
    # --endwndw       2022-11-02_00:00:00.0 \
    # --plot          rs \

python ~/scripts/seafloor-geod/apg_read.py \
    --infile        ./TAN2301_GNS22-PO_140342.apg \
    --snapg         140342 \
    --version       CSAC2013 \
    --tempsmth      5001 \
    --decimate      1 \
    --clkstart      2022-10-11_05:08:00 \
    --bininterval   1d \
    --gpssynctime   2023:007:06:02:00 \
    --synctickcount 0x01C561104A \
    --fmttime       d \
    --outfile       ./out/GNS22-PO.csv \
    | tee ./out/GNS22-PO.log
    # --beginwndw     2022-11-01_00:00:00.0 \
    # --endwndw       2022-11-02_00:00:00.0 \
    # --plot          rs \
