# Created : Thanos Tsikerdekis | Apr 2025 | Download from MARS
# Parts of code from J.Douros and B.Mijling scripts were used and modified.
# Designed for M7 pressure (pl) levels

####################
### LOAD MODULES ###
####################
import os, sys
import subprocess
import pandas as pd
from ecmwfapi import ECMWFService
from datetime import datetime, timedelta, date
from dateutil import parser
from tenacity import *
from random   import randint

################################
### FUNCTION retrieve_global ###
################################
def retrieve_global(expname,expclass,day):
    
    ### SET PATHS
    TEMP_DIR  = '/tmp/'
    OUT_DIR   = path_data + expname
    DATA_NAME = 'CAMS_' + expname + '_forecast00to21by03_0.7x0.7_pl'

    ### CREATE TEMPORARY DIRECTORY IF IT DOES NOT EXIST
    daystrip = day.replace('-','')
    mytempdir = TEMP_DIR + daystrip
    if not os.path.exists(mytempdir): os.makedirs(mytempdir)

    # CREATE OUTPUT DIRECTORY IF IT DOES NOT EXIST
    if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)
    myoutput = OUT_DIR + "/" + DATA_NAME + "_" + daystrip + ".nc"

#    # SKIP IF OUTPUT FILE ALREADY EXISTS
#    if os.path.isfile(myoutput):
#        print(f"File {myoutput} already exists, skipping...")
#        return #Skip downloading if file exists

    # IF TEMPORARY FILE EXIST REMOVE IT
    mytempfn = TEMP_DIR + "Temp_" + DATA_NAME + "_" + day.replace('-','') + "_day_3D.nc"
    if os.path.isfile(mytempfn):
        command = 'rm -f ' + mytempfn
        subprocess.run(command, shell=True)

#    @retry(stop=stop_after_attempt(1))
#    def retryfc():
#        if True:
#            print("Trying to download:",mytempfn)
#            server = ECMWFService("mars")
#            server.execute({
#                "class": expclass,
#                "date": day,
#                "expver": expname,
#                "levelist": "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000",
#                "levtype": "pl",
#                # Geopotential = 129.128 # pl
#                # Temperature  = 130.128 # pl
#                "param": "1.212/2.212/3.212/4.212/5.212/6.212/7.212/8.212/9.212/10.212/11.212/12.212/13.212/14.212/15.212/16.212/17.212/18.212/19.212/20.212/21.212/22.212/23.212/24.212/25.212/26.212/27.212/28.212/29.212/30.212/31.212/32.212/33.212/34.212/35.212/36.212/37.212/38.212/129.128/130.128",
#                "step": "0/3/6/9/12/15/18/21",
#                "stream": "oper",
#                "time": "00",
#                "type": "fc",
#                "format": "netcdf",
#                "grid": "0.7/0.7",
#                },
#                mytempfn)
#            print("Apparently succeded!")
#
#        else:
#            #except Exception as e:
#            print("The following error happened while trying to download:", str(e))
#            # If the following exception happen, do not retry
#            if ('No available data matches request' or \
#                    'Bad magic number') in str(e):
#                    print("This is a special kind of error, stop retrying...")
#                    return
#    retryfc()
#
#    # COMPRESS
#    command = 'ncpdq -O -4 -L 1 ' + mytempfn + ' ' + myoutput
#    subprocess.run(command, shell=True)

    # RENAME M7 VARIABLES
    command = (
    "ncrename "
    "-v p212001,AS_N "
    "-v p212002,SO4_AS "
    "-v p212003,BC_AS "
    "-v p212004,POM_AS "
    "-v p212005,SS_AS "
    "-v p212006,DU_AS "
    "-v p212007,SOA_NS "
    "-v p212008,SOA_KS "
    "-v p212009,SOA_AS "
    "-v p212010,SOA_CS "
    "-v p212011,SOA_KI "
    "-v p212012,H2OPART "
    "-v p212013,KI_N "
    "-v p212014,BC_KI "
    "-v p212015,POM_KI "
    "-v p212016,AI_N "
    "-v p212017,DU_AI "
    "-v p212018,KS_N "
    "-v p212019,SO4_KS "
    "-v p212020,BC_KS "
    "-v p212021,POM_KS "
    "-v p212022,CI_N "
    "-v p212023,DU_CI "
    "-v p212024,CS_N "
    "-v p212025,SO4_CS "
    "-v p212026,BC_CS "
    "-v p212027,POM_CS "
    "-v p212028,SS_CS "
    "-v p212029,DU_CS "
    "-v p212030,NS_N "
    "-v p212031,SO4_NS "
    "-v p212032,ELVOC "
    "-v p212033,ISVOC "
    "-v p212034,MSA "
    "-v p212035,NH4 "
    "-v p212036,NO3_A "
    "-v p212037,CDNC "
    "-v p212038,ICNC "
    + myoutput
    )
    subprocess.run(command, shell=True)

    # MIXING RATIO
    command = "ncap2 -O -s 'MR_DU=DU_AI+DU_CI+DU_AS+DU_CS;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_SS=SS_AS+SS_CS;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_POM=POM_KI+POM_KS+POM_AS+POM_CS;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_BC=BC_KI+BC_KS+BC_AS+BC_CS;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_SO4=SO4_NS+SO4_KS+SO4_AS+SO4_CS;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_NI=NO3_A;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_AM=NH4;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_AERO=MR_DU+MR_SS+MR_POM+MR_BC+MR_SO4+MR_NI+MR_AM;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    
    # TOTAL COLUMN BURDEN
    species_list = ["DU", "SS", "POM", "BC", "SO4", "NI", "AM", "AERO"]
    for sp in species_list:
      command = (
          f"ncap2 -O -s '"
          f"g=9.80665f; "
          f"dp[level]={{50.f,100.f,150.f,200.f,250.f,650.f,1000.f,1500.f,2000.f,2500.f,4000.f,5000.f,5000.f,5000.f,7500.f,10000.f,15000.f,17500.f,11250.f,7500.f,5075.f}}; "
          f"BU_{sp}[$time,$latitude,$longitude]=0.0f; "
          f"for (l=0; l<21; l++) {{ BU_{sp}(:,:,:) += MR_{sp}(:,l,:,:) * dp(l)/g; }}'"
          f" {myoutput} -o {myoutput}")
      subprocess.run(command, shell=True)

    # CLEANUP TEMP FILES
    command = "rm -f " + mytempfn
    subprocess.run(command, shell=True)

    # SHARE THE GOOD NEWS
    print("File output at:", myoutput)

############
### MAIN ###
############
if __name__ == "__main__":

    if "-h" in sys.argv or "--help" in sys.argv:
        print("\nDescription : Download daily files from a specific experiment using MARS")
        print("Usage       : python script.py EXPNAME DATESTART DATEEND")
        print("Example1    : python script.py b2r3 20250101")
        print("Example2    : python script.py b2r3 20250101 20250103\n")
        sys.exit()

    # If system argument is one: Download model data some days ago
    # The idea is to give some extra time for all the models to upload files
    if len(sys.argv)!=4:
        print("\nToo few or too many arguments. The script needs 2 or 3 arguments.")
        print("\nDescription : Download daily files from a specific experiment using MARS")
        print("Usage       : python script.py EXPNAME DATESTART DATEEND")
        print("Example1    : python script.py b2r3 20250101")
        print("Example2    : python script.py b2r3 20250101 20250103\n")
        sys.exit()

    # If system arguments are four: Download data file for days between this dates + define output (for COMPASS)
    # Each day count as on MARS request
    # Provide them as: python download_crontab_CAMSr.py 20220627 20220630 path_data
    elif len(sys.argv)==4:
        expname      = sys.argv[1]
        startDate    = datetime.strptime(sys.argv[2],'%Y%m%d')
        endDate      = datetime.strptime(sys.argv[2],'%Y%m%d')
        sequenceDate = pd.date_range(startDate, endDate, freq='D')
        dateCodes    = sequenceDate.strftime('%Y-%m-%d')
        path_data    = sys.argv[3]

    ### Loop through dates and submit download requests for control and analysis
    for dateCode in dateCodes:
        retrieve_global(expname,"nl",dateCode) # osuite

