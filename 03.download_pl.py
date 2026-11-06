#!/usr/bin/env python3
"""PL 30-day GRIB downloader; keeps legacy daily COMPASS NetCDF names."""
import sys
from download_batch import retrieve_batches

PL_LEVELS = "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000"


def main(args):
    if len(args)!=6:
        raise SystemExit("Usage: 03.download_pl.py EXPNAME EXPCLASS DATESTART DATEEND PATH_DATA PARAMS")
    expname,expclass,start,end,path_data,params=args
    retrieve_batches(expname,expclass,start,end,path_data,params,
                     tag="pl",levtype="pl",levels=PL_LEVELS)


if __name__=="__main__":
    main(sys.argv[1:])
