"""Shared MARS GRIB batch retrieval and safe daily NetCDF conversion for COMPASS.

Requires ECMWFService, CDO >= 2.x, NCO ncrename, and ncdump on the cluster.
Daily outputs retain the existing COMPASS names and initialisation-day groups.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

from ecmwfapi import ECMWFService

CDO = os.environ.get("COMPASS_CDO", "/usr/local/apps/cdo/2.5.1/bin/cdo")
STEPS_PER_DAY = 8
MAX_BATCH_DAYS = 30


def run(*args):
    cmd = [str(x) for x in args]
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()


def days_inclusive(start, end):
    a = datetime.strptime(start, "%Y%m%d").date()
    b = datetime.strptime(end, "%Y%m%d").date()
    if b < a:
        raise ValueError("DATEEND must be >= DATESTART")
    return [a + timedelta(days=i) for i in range((b-a).days+1)]


def groups(dates, max_days=MAX_BATCH_DAYS):
    """Only contiguous missing days can belong to the same MARS request."""
    pending = []
    for day in dates:
        if pending and (day != pending[-1]+timedelta(days=1) or len(pending) == max_days):
            yield pending
            pending = []
        pending.append(day)
    if pending:
        yield pending


def canonical(code):
    num, table = code.split(".")
    return f"p{int(table):03d}{int(num):03d}"


def normalized_code(code):
    val = code.strip()
    if "." in val:
        n,t=val.split(".")
        return int(t)*1000+int(n)
    return int(val)


def check_existing(path, params=None):
    """A nonempty existing NetCDF is preserved; validate contents before reuse."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        ntime = int(run(CDO, "-s", "ntime", path).split()[-1])
        if ntime != STEPS_PER_DAY:
            return False
        if params:
            header=nc_header(path)
            for param in params.split("/"):
                if not re.search(r"\b"+canonical(param)+r"\s*\(",header):
                    return False
        return True
    except (OSError, ValueError, subprocess.CalledProcessError):
        return False


def nc_header(path):
    return run("ncdump", "-h", path)


def rename_coordinates(path):
    """CDO often writes lon/lat/plev instead of MARS longitude/latitude/level."""
    head = nc_header(path)
    for old,new in (("lon","longitude"),("lat","latitude"),
                    ("plev","level"),("lev","level")):
        if re.search(r"\b"+new+r"\s*=", head):
            continue
        if not re.search(r"\b"+old+r"\s*=", head):
            continue
        flags=["ncrename","-O","-d",f"{old},{new}"]
        if re.search(r"(?:^|\s)\w+\s+"+old+r"\s*\(",head):
            flags.extend(["-v",f"{old},{new}"])
        run(*flags,path)
        head=nc_header(path)


def normalize_variables(grib, nc, params):
    """Map actual GRIB parameter codes to canonical COMPASS variable names."""
    wanted={normalized_code(x): canonical(x) for x in params.split("/")}
    actual=run(CDO,"-s","showparam",grib).split()
    source_names=run(CDO,"-s","showname",grib).split()
    names=run(CDO,"-s","showname",nc).split()
    print("GRIB parameters:",actual,"GRIB names:",source_names,
          "NetCDF variables:",names,flush=True)
    if source_names != names:
        raise RuntimeError("CDO changed variable ordering or names in GRIB-to-NetCDF "
                           f"conversion. Refusing ambiguous renaming: {source_names} -> {names}")
    if len(actual) != len(names):
        raise RuntimeError("Parameter/name count mismatch: cannot safely rename " + str((actual,names)))
    found={}
    for param,name in zip(actual,names):
        code=normalized_code(param)
        if code not in wanted or code in found:
            raise RuntimeError(f"Unexpected or repeated GRIB parameter {param} in {grib}")
        found[code]=name
    if set(found) != set(wanted):
        raise RuntimeError(f"Missing parameters: {set(wanted)-set(found)}")
    # Rename in two passes so a source name cannot shadow a destination name.
    for i,(code,name) in enumerate(found.items()):
        if name != wanted[code]:
            run("ncrename","-O","-v",f"{name},_compass_tmp_{i}",nc)
            found[code]=f"_compass_tmp_{i}"
    for code,name in found.items():
        if name != wanted[code]:
            run("ncrename","-O","-v",f"{name},{wanted[code]}",nc)
    rename_coordinates(nc)


def verify_daily_times(grib, day, steps):
    """Verify the eight GRIB valid times match the requested forecast init day.

    In particular step +24 belongs to the requested day's output even though
    its valid time falls on the next calendar day.
    """
    values=run(CDO,"-s","showtimestamp",grib).split()
    if len(values)!=STEPS_PER_DAY:
        raise RuntimeError(f"Expected eight valid timestamps for {day}: {values}")
    hours=[int(x) for x in steps.split("/")]
    if len(hours)!=STEPS_PER_DAY:
        raise RuntimeError(f"Expected eight forecast steps, got {steps}")
    init=datetime.combine(day,datetime.min.time())
    expected=[init+timedelta(hours=h) for h in hours]
    try:
        actual=[datetime.fromisoformat(x.replace("Z","+00:00")).replace(tzinfo=None)
                for x in values]
    except ValueError as err:
        raise RuntimeError(f"Unrecognized CDO timestamp output: {values}") from err
    if actual!=expected:
        raise RuntimeError(f"GRIB times for initialization day {day} are wrong: "
                           f"got {actual}, expected {expected}")


def convert_one(grib, target, params):
    """Never publish partially converted output."""
    tmp=target.with_name(target.name+f".partial.{os.getpid()}")
    try:
        run(CDO,"-s","-f","nc4","-z","zip_1","copy",grib,tmp)
        normalize_variables(grib,tmp,params)
        ntime=int(run(CDO,"-s","ntime",tmp).split()[-1])
        if ntime != STEPS_PER_DAY:
            raise RuntimeError(f"Expected {STEPS_PER_DAY} timesteps but found {ntime} in {tmp}")
        header=nc_header(tmp)
        if not re.search(r"\blongitude\s*=",header) or not re.search(r"\blatitude\s*=",header):
            raise RuntimeError(f"Missing COMPASS longitude/latitude dimensions in {tmp}")
        if not target.exists() or not check_existing(target,params):
            os.replace(tmp,target)
            print("Completed:",target,flush=True)
        else:
            print("Another worker already finished:",target,flush=True)
    finally:
        tmp.unlink(missing_ok=True)


def retrieve_batches(expname,expclass,start,end,path_data,params,*,tag,levtype,
                     levels=None,steps="0/3/6/9/12/15/18/21"):
    out_dir=Path(path_data)/expname
    out_dir.mkdir(parents=True,exist_ok=True)
    prefix=(f"CAMS_{expname}_forecast{int(steps.split('/')[0]):02d}"
            f"to{int(steps.split('/')[-1]):02d}by03_0.7x0.7_{tag}")
    dates=days_inclusive(start,end)
    missing=[d for d in dates if not check_existing(out_dir/f"{prefix}_{d:%Y%m%d}.nc",params)]
    print(f"{tag}: {len(dates)-len(missing)} existing valid days, "
          f"{len(missing)} missing/invalid days",flush=True)
    if not missing:
        return
    server=ECMWFService("mars")
    for chunk in groups(missing):
        first,last=chunk[0],chunk[-1]
        print(f"MARS {tag}: requesting {first} through {last} ({len(chunk)} days)",flush=True)
        # Temporary files are unique to the job and reside on /scratch, not /tmp.
        with tempfile.TemporaryDirectory(prefix=f".COMPASS_{tag}_",dir=out_dir) as scratch:
            work=Path(scratch)
            batch_grib=work/"batch.grib"
            request={
                "class":expclass,"date":f"{first:%Y-%m-%d}/to/{last:%Y-%m-%d}",
                "expver":expname,"levtype":levtype,"param":params,
                "step":steps,"stream":"oper","time":"00","type":"fc",
                "grid":"0.7/0.7"
            }
            if levels:
                request["levelist"]=levels
            print("MARS request:",request,flush=True)
            server.execute(request,str(batch_grib))
            if not batch_grib.is_file() or not batch_grib.stat().st_size:
                raise RuntimeError("MARS returned an empty GRIB: "+str(batch_grib))
            ntime=int(run(CDO,"-s","ntime",batch_grib).split()[-1])
            if ntime != len(chunk)*STEPS_PER_DAY:
                raise RuntimeError(f"MARS GRIB has {ntime} times; expected {len(chunk)*STEPS_PER_DAY}. "
                                   "Keeping the existing daily NetCDF untouched.")
            # Eight initialisation-time steps per day; do NOT use splitdate,
            # since precipitation includes step 24 (valid on the following day).
            split_prefix=work/"day_"
            run(CDO,"-s",f"splitsel,{STEPS_PER_DAY}",batch_grib,split_prefix)
            pieces=sorted(work.glob("day_*"))
            if len(pieces) != len(chunk):
                raise RuntimeError(f"CDO generated {len(pieces)} daily GRIBs for {len(chunk)} days")
            for day,grib in zip(chunk,pieces):
                target=out_dir/f"{prefix}_{day:%Y%m%d}.nc"
                if check_existing(target,params):
                    print("Exists; skip:",target,flush=True)
                    continue
                day_nt=int(run(CDO,"-s","ntime",grib).split()[-1])
                if day_nt != STEPS_PER_DAY:
                    raise RuntimeError(f"Split file for {day} has {day_nt} timesteps")
                verify_daily_times(grib,day,steps)
                convert_one(grib,target,params)
