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
    DATA_NAME = 'CAMS_' + expname + '_forecast00to21by03_0.7x0.7_sfc'

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
    mytempfn = TEMP_DIR + "Temp_" + DATA_NAME + "_" + day.replace('-','') + "_day_2D.nc"
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
                "levtype": "sfc",
                "param": "104.215/109.215/140.215/207.210/215.210",
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
    if len(sys.argv)!=6:
        print("\nToo few or too many arguments. The script needs 2 or 3 arguments.")
        print("\nDescription : Download daily files from a specific experiment using MARS")
        print("Usage       : python script.py EXPNAME DATESTART DATEEND")
        print("Example1    : python script.py b2r3 20250101")
        print("Example2    : python script.py b2r3 20250101 20250103\n")
        sys.exit()

    # If system arguments are four: Download data file for days between this dates + define output (for COMPASS)
    # Each day count as on MARS request
    # Provide them as: python download_crontab_CAMSr.py 20220627 20220630 path_data
    elif len(sys.argv)==6:
        expname      = sys.argv[1]
        expclass     = sys.argv[2]
        startDate    = datetime.strptime(sys.argv[3],'%Y%m%d')
        endDate      = datetime.strptime(sys.argv[4],'%Y%m%d')
        sequenceDate = pd.date_range(startDate, endDate, freq='D')
        dateCodes    = sequenceDate.strftime('%Y-%m-%d')
        path_data    = sys.argv[5]

    ### Loop through dates and submit download requests for experiment
    for dateCode in dateCodes:
        retrieve_global(expname,expclass,dateCode) 

