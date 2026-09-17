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

    print("Requested GRIBs :",requested)
    print("Current variables:",current_names)
    print("Expected names   :",expected_names)

    if len(current_names) != len(expected_names):
        raise RuntimeError(
            f"Cannot safely rename variables. Requested {len(expected_names)} fields "
            f"but NetCDF contains {len(current_names)} data variables: {current_names}"
        )

    for current,expected in zip(current_names,expected_names):
        if current == expected:
            print(f"Variable {current} already correctly named.")
            continue
        print(f"Renaming {current} -> {expected}")
        subprocess.run(["ncrename","-O","-v",f"{current},{expected}",filename],check=True)


def retrieve_global(expname,expclass,day,path_data,levels):
    TEMP_DIR = "/tmp/"
    OUT_DIR = path_data + expname
    daystrip = day.replace("-","")

    if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)

    # T and q are stored one model level per file. This avoids ambiguity when
    # different RH heights are requested in separate COMPASS runs.
    params_ml = "130.128/133.128"
    model_levels = sorted({int(x) for x in levels.split("/")})

    server = ECMWFService("mars")

    for model_level in model_levels:
        output_ml = OUT_DIR + f"/CAMS_{expname}_forecast00to21by03_0.7x0.7_ml_{model_level}_{daystrip}.nc"
        temp_ml = TEMP_DIR + f"Temp_CAMS_{expname}_ml_{model_level}_{daystrip}.nc"

        if os.path.isfile(output_ml):
            print(f"File {output_ml} already exists, skipping...")
            continue

        if os.path.isfile(temp_ml): os.remove(temp_ml)

        @retry(stop=stop_after_attempt(1))
        def retry_ml():
            print("Trying to download model-level T/q:",temp_ml)
            print("Model level:",model_level)

            server.execute({
                "class": expclass,
                "date": day,
                "expver": expname,
                "levelist": str(model_level),
                "levtype": "ml",
                "param": params_ml,
                "step": "0/3/6/9/12/15/18/21",
                "stream": "oper",
                "time": "00",
                "type": "fc",
                "format": "netcdf",
                "grid": "0.7/0.7",
            },temp_ml)

        retry_ml()
        normalize_variable_names(temp_ml,params_ml)
        subprocess.run(["ncpdq","-O","-4","-L","1",temp_ml,output_ml],check=True)

        if os.path.isfile(temp_ml): os.remove(temp_ml)
        print("Model-level output:",output_ml)

    # ln(surface pressure), GRIB 152, is shared by all requested RH levels.
    params_lnsp = "152.128"
    output_lnsp = OUT_DIR + f"/CAMS_{expname}_forecast00to21by03_0.7x0.7_lnsp_{daystrip}.nc"
    temp_lnsp = TEMP_DIR + f"Temp_CAMS_{expname}_lnsp_{daystrip}.nc"

    if not os.path.isfile(output_lnsp):
        if os.path.isfile(temp_lnsp): os.remove(temp_lnsp)

        @retry(stop=stop_after_attempt(1))
        def retry_lnsp():
            print("Trying to download lnsp:",temp_lnsp)

            server.execute({
                "class": expclass,
                "date": day,
                "expver": expname,
                "levelist": "1",
                "levtype": "ml",
                "param": params_lnsp,
                "step": "0/3/6/9/12/15/18/21",
                "stream": "oper",
                "time": "00",
                "type": "fc",
                "format": "netcdf",
                "grid": "0.7/0.7",
            },temp_lnsp)

        retry_lnsp()
        normalize_variable_names(temp_lnsp,params_lnsp)
        subprocess.run(["ncpdq","-O","-4","-L","1",temp_lnsp,output_lnsp],check=True)

        if os.path.isfile(temp_lnsp): os.remove(temp_lnsp)
        print("lnsp output:",output_lnsp)
    else:
        print(f"File {output_lnsp} already exists, skipping...")


if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        print("\nUsage: python 03.download_ml.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA LEVELS\n")
        sys.exit()

    if len(sys.argv) != 7:
        print("\nIncorrect number of arguments.")
        print("Usage: python 03.download_ml.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA LEVELS\n")
        sys.exit(1)

    expname = sys.argv[1]
    expclass = sys.argv[2]
    startDate = datetime.strptime(sys.argv[3],"%Y%m%d")
    endDate = datetime.strptime(sys.argv[4],"%Y%m%d")
    path_data = sys.argv[5]
    levels = sys.argv[6]

    sequenceDate = pd.date_range(startDate,endDate,freq="D")
    dateCodes = sequenceDate.strftime("%Y-%m-%d")

    for dateCode in dateCodes:
        retrieve_global(expname,expclass,dateCode,path_data,levels)
