#!/anaconda3/bin/python -u

import os
import sys
import subprocess
import pandas as pd
from ecmwfapi import ECMWFService
from datetime import datetime
from tenacity import retry, stop_after_attempt

CDO = "/usr/local/apps/cdo/2.5.1/bin/cdo"

def grib_to_canonical_ncname(grib):
    number, table = grib.split(".")
    return f"p{int(table):03d}{int(number):03d}"

def normalize_variable_names(filename, params):

    requested = params.split("/")
    expected_names = [grib_to_canonical_ncname(x) for x in requested]

    result = subprocess.run([CDO,"-s","showname",filename],capture_output=True,text=True,check=True)
    current_names = result.stdout.split()

    print("Requested GRIBs :", requested)
    print("Current variables:", current_names)
    print("Expected names   :", expected_names)

    if len(current_names) != len(expected_names):
        raise RuntimeError(
            f"Cannot safely rename variables. Requested {len(expected_names)} fields "
            f"but NetCDF contains {len(current_names)} data variables: {current_names}"
        )

    for current, expected in zip(current_names,expected_names):

        if current == expected:
            print(f"Variable {current} already correctly named.")
            continue

        print(f"Renaming {current} -> {expected}")
        subprocess.run(["ncrename","-O","-v",f"{current},{expected}",filename],check=True)

def retrieve_global(expname, expclass, day, params):

    TEMP_DIR = "/tmp/"
    OUT_DIR = path_data + expname
    DATA_NAME = "CAMS_" + expname + "_forecast00to21by03_0.7x0.7_sfc"

    daystrip = day.replace("-","")
    mytempdir = TEMP_DIR + daystrip

    if not os.path.exists(mytempdir): os.makedirs(mytempdir)
    if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)

    myoutput = OUT_DIR + "/" + DATA_NAME + "_" + daystrip + ".nc"

    if os.path.isfile(myoutput):
        print(f"File {myoutput} already exists, skipping...")
        return

    mytempfn = TEMP_DIR + "Temp_" + DATA_NAME + "_" + daystrip + "_day_2D.nc"
    if os.path.isfile(mytempfn): os.remove(mytempfn)

    @retry(stop=stop_after_attempt(1))
    def retryfc():

        print("Trying to download:",mytempfn)
        print("Parameters:",params)

        server = ECMWFService("mars")

        server.execute({
            "class": expclass,
            "date": day,
            "expver": expname,
            "levtype": "sfc",
            "param": params,
            "step": "0/3/6/9/12/15/18/21",
            "stream": "oper",
            "time": "00",
            "type": "fc",
            "format": "netcdf",
            "grid": "0.7/0.7",
        },mytempfn)

        print("Apparently succeeded!")

    retryfc()

    normalize_variable_names(mytempfn,params)

    subprocess.run(["ncpdq","-O","-4","-L","1",mytempfn,myoutput],check=True)

    if os.path.isfile(mytempfn): os.remove(mytempfn)

    print("File output at:",myoutput)

if __name__ == "__main__":

    if "-h" in sys.argv or "--help" in sys.argv:
        print("\nUsage: python script.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA PARAMS\n")
        sys.exit()

    if len(sys.argv) != 7:
        print("\nIncorrect number of arguments.")
        print("Usage: python script.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA PARAMS\n")
        sys.exit(1)

    expname = sys.argv[1]
    expclass = sys.argv[2]
    startDate = datetime.strptime(sys.argv[3],"%Y%m%d")
    endDate = datetime.strptime(sys.argv[4],"%Y%m%d")
    path_data = sys.argv[5]
    params = sys.argv[6]

    sequenceDate = pd.date_range(startDate,endDate,freq="D")
    dateCodes = sequenceDate.strftime("%Y-%m-%d")

    for dateCode in dateCodes:
        retrieve_global(expname,expclass,dateCode,params)
