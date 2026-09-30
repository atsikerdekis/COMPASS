#!/anaconda3/bin/python -u

import os
import sys
import subprocess
from ecmwfapi import ECMWFService
from datetime import datetime, timedelta
from tenacity import retry, stop_after_attempt

CDO = "/usr/local/apps/cdo/2.5.1/bin/cdo"


def grib_to_canonical_ncname(grib):
    number, table = grib.split(".")
    return f"p{int(table):03d}{int(number):03d}"


def normalize_variable_names(filename, params):
    requested = params.split("/")
    expected_names = [grib_to_canonical_ncname(x) for x in requested]

    result = subprocess.run([CDO, "-s", "showname", filename], capture_output=True, text=True, check=True)
    current_names = result.stdout.split()

    print("Requested GRIBs :", requested)
    print("Current variables:", current_names)
    print("Expected names   :", expected_names)

    if len(current_names) != len(expected_names):
        raise RuntimeError(
            f"Cannot safely rename variables. Requested {len(expected_names)} fields "
            f"but NetCDF contains {len(current_names)} data variables: {current_names}"
        )

    for current, expected in zip(current_names, expected_names):
        if current == expected:
            print(f"Variable {current} already correctly named.")
            continue
        print(f"Renaming {current} -> {expected}")
        subprocess.run(["ncrename", "-O", "-v", f"{current},{expected}", filename], check=True)


def mars_to_netcdf(server, request, temp_file, out_file, params):
    if os.path.isfile(out_file):
        print(f"File {out_file} already exists, skipping...")
        return
    if os.path.isfile(temp_file):
        os.remove(temp_file)

    @retry(stop=stop_after_attempt(1))
    def _run():
        print("Trying to download:", temp_file)
        print("Request:", request)
        server.execute(request, temp_file)

    _run()
    normalize_variable_names(temp_file, params)
    subprocess.run(["ncpdq", "-O", "-4", "-L", "1", temp_file, out_file], check=True)
    if os.path.isfile(temp_file):
        os.remove(temp_file)
    print("Output:", out_file)


def retrieve_global(expname, expclass, day, path_data, levels, params_ml="130.128/133.128", tag="rh"):
    temp_dir = "/tmp/"
    out_dir = os.path.join(path_data, expname)
    daystrip = day.replace("-", "")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    server = ECMWFService("mars")

    output_ml = os.path.join(out_dir, f"CAMS_{expname}_forecast00to21by03_0.7x0.7_ml_{tag}_{daystrip}.nc")
    temp_ml = os.path.join(temp_dir, f"Temp_CAMS_{expname}_ml_{tag}_{daystrip}.nc")

    req_ml = {
        "class": expclass,
        "date": day,
        "expver": expname,
        "levelist": levels,
        "levtype": "ml",
        "param": params_ml,
        "step": "0/3/6/9/12/15/18/21",
        "stream": "oper",
        "time": "00",
        "type": "fc",
        "format": "netcdf",
        "grid": "0.7/0.7",
    }
    mars_to_netcdf(server, req_ml, temp_ml, output_ml, params_ml)

    if tag == "rh":
        params_lnsp = "152.128"
        output_lnsp = os.path.join(out_dir, f"CAMS_{expname}_forecast00to21by03_0.7x0.7_lnsp_{daystrip}.nc")
        temp_lnsp = os.path.join(temp_dir, f"Temp_CAMS_{expname}_lnsp_{daystrip}.nc")
        req_lnsp = {
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
        }
        mars_to_netcdf(server, req_lnsp, temp_lnsp, output_lnsp, params_lnsp)


if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        print("\nUsage: python 03.download_ml.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA LEVELS [PARAMS] [TAG]\n")
        sys.exit()

    if len(sys.argv) < 7 or len(sys.argv) > 9:
        print("\nIncorrect number of arguments.")
        print("Usage: python 03.download_ml.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA LEVELS [PARAMS] [TAG]\n")
        sys.exit(1)

    expname = sys.argv[1]
    expclass = sys.argv[2]
    start_date = datetime.strptime(sys.argv[3], "%Y%m%d")
    end_date = datetime.strptime(sys.argv[4], "%Y%m%d")
    path_data = sys.argv[5]
    levels = sys.argv[6]
    params_ml = sys.argv[7] if len(sys.argv) >= 8 else "130.128/133.128"
    tag = sys.argv[8] if len(sys.argv) >= 9 else "rh"

    date = start_date
    while date <= end_date:
        retrieve_global(expname, expclass, date.strftime("%Y-%m-%d"), path_data, levels, params_ml=params_ml, tag=tag)
        date += timedelta(days=1)
