#!/usr/bin/env python3
"""ML 30-day GRIB downloader, including standalone lnsp for RH."""
import sys
from download_batch import retrieve_batches


def main(args):
    if len(args) not in (6,7,8):
        raise SystemExit("Usage: 03.download_ml.py EXPNAME EXPCLASS DATESTART DATEEND "
                         "PATH_DATA LEVELS [PARAMS] [TAG]")
    expname,expclass,start,end,path_data,levels=args[:6]
    params=args[6] if len(args)>=7 else "130.128/133.128"
    tag=args[7] if len(args)>=8 else "rh"
    retrieve_batches(expname,expclass,start,end,path_data,params,
                     tag=f"ml_{tag}",levtype="ml",levels=levels)
    if tag=="rh":
        retrieve_batches(expname,expclass,start,end,path_data,"152.128",
                         tag="lnsp",levtype="ml",levels="1")


if __name__=="__main__":
    main(sys.argv[1:])
