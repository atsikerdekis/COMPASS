#!/usr/bin/env python3
"""SFC 30-day GRIB downloader; keeps legacy daily COMPASS NetCDF names."""
import sys
from datetime import datetime
from download_batch import retrieve_batches


def main(args):
    if len(args) not in (6, 7, 8):
        raise SystemExit("Usage: 03.download_sfc.py EXPNAME EXPCLASS DATESTART DATEEND "
                         "PATH_DATA PARAMS [STEPS] [TAG]")
    expname,expclass,start,end,path_data,params=args[:6]
    steps=args[6] if len(args)>=7 else "0/3/6/9/12/15/18/21"
    tag=args[7] if len(args)>=8 else "sfc"
    datetime.strptime(start,"%Y%m%d")
    datetime.strptime(end,"%Y%m%d")
    retrieve_batches(expname,expclass,start,end,path_data,params,
                     tag=tag,levtype="sfc",steps=steps)


if __name__=="__main__":
    main(sys.argv[1:])
