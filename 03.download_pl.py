#!/anaconda3/bin/python -u

###!/anaconda3/envs/cams2_83/bin/python

# Created : Thanos Tsikerdekis | Apr 2025 | Download from MARS
# Parts of code from J.Douros and B.Mijling scripts were used and modified.

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

    # SKIP IF OUTPUT FILE ALREADY EXISTS
    if os.path.isfile(myoutput):
        print(f"File {myoutput} already exists, skipping...")
        return #Skip downloading if file exists

    # IF TEMPORARY FILE EXIST REMOVE IT
    mytempfn = TEMP_DIR + "Temp_" + DATA_NAME + "_" + day.replace('-','') + "_day_3D.nc"
    if os.path.isfile(mytempfn):
        command = 'rm -f ' + mytempfn
        subprocess.run(command, shell=True)

    @retry(stop=stop_after_attempt(1))
    def retryfc():
        if True:
            print("Trying to download:",mytempfn)
            server = ECMWFService("mars")
            server.execute({
                "class": expclass,
                "date": day,
                "expver": expname,
                "levelist": "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000",
                "levtype": "pl",
                # Geopotential = 129.128 # pl
                # Temperature  = 130.128 # pl
                "param": "1.210/2.210/3.210/4.210/5.210/6.210/7.210/8.210/9.210/10.210/11.210/247.210/248.210/249.210/210252/210253/129.128/130.128",
                "step": "0/3/6/9/12/15/18/21",
                "stream": "oper",
                "time": "00",
                "type": "fc",
                "format": "netcdf",
                "grid": "0.7/0.7",
                },
                mytempfn)
            print("Apparently succeded!")

        else:
            #except Exception as e:
            print("The following error happened while trying to download:", str(e))
            # If the following exception happen, do not retry
            if ('No available data matches request' or \
                    'Bad magic number') in str(e):
                    print("This is a special kind of error, stop retrying...")
                    return
    retryfc()

    # COMPRESS
    command = 'ncpdq -O -4 -L 1 ' + mytempfn + ' ' + myoutput
    subprocess.run(command, shell=True)

    # MIXING RATIO
    command = "ncap2 -O -s 'MR_DU=aermr04+aermr05+aermr06;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_SS=aermr01+aermr02+aermr03;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_POM=aermr07+aermr08+aermr19+aermr20;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_BC=aermr09+aermr10;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_SO4=aermr11;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_NI=aermr16+aermr17;' " + myoutput + " -o " + myoutput
    subprocess.run(command, shell=True)
    command = "ncap2 -O -s 'MR_AM=aermr18;' " + myoutput + " -o " + myoutput
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

