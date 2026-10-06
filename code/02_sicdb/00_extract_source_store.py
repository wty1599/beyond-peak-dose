"""Extract authorized local SICdb data to a bounded SQLite derived store.

Raw files are never modified. The store is restricted patient-level working data,
not a public deliverable. No coefficients, imputation, predictions or inference.
"""
from pathlib import Path
from collections import Counter
import argparse
import json
import sqlite3
import sys
import time
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'00_setup'))
from paths import release_path, private_output
DATA = release_path('sicdb_raw')
OUT = private_output('sicdb_source_store')
WORK = OUT
NA = ["", "NULL", "null", "\\N", "NA", "NaN"]
CON = sqlite3.connect(WORK / "sicdb_selected.sqlite")
CON.execute("PRAGMA cache_size=-65536")
CON.execute("PRAGMA temp_store=FILE")
CHUNK = 100000
DRUGS = {1502, 1550, 1562, 1593, 1618, 1559, 1560, 2046}
SIGNALS = {703,706,707,708,724,725,722,723,730,731,732,2022,3071,
           711,712,713,715,716,717,718,719,727,2019,2020,2021,
           2023,2024,2025,2280,2281,2278,2282,2283,2284,2285,3040,
           3139,3123,710}
LABS = {454,465,657,314,333,367,664,689}

def log(msg):
    print(time.strftime("%H:%M:%S"),msg,flush=True)

def read(name, **kw):
    p=DATA/(name+".csv.gz")
    assert p.is_file() and not p.is_symlink()
    return pd.read_csv(p,na_values=NA,**kw)

def export_json(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")

def already(table):
    return CON.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",(table,)).fetchone() is not None

def prep():
    if already("medication"):
        raise RuntimeError("Preparation output already exists; no automatic overwrite")
    c=read("cases")
    c.to_sql("cases",CON,index=False)
    d=read("d_references")
    d.to_sql("dictionary",CON,index=False)
    for table in ["data_ref","data_range","unitlog"]:
        a=read(table)
        a=a[a.CaseID.isin(c.CaseID)]
        a.to_sql(table,CON,index=False)
    valid=c[(c.AgeOnAdmission>=18)&(c.TimeOfStay>c.ICUOffset)&(c.WeightOnAdmission>0)]
    ids=set(valid.CaseID.astype(int))
    parts=[]
    medrows=0
    for chunk in read("medication",chunksize=CHUNK):
        medrows+=len(chunk)
        keep=chunk[chunk.CaseID.isin(ids)&chunk.DrugID.isin(DRUGS)]
        if len(keep): parts.append(keep)
    m=pd.concat(parts,ignore_index=True)
    m.to_sql("medication",CON,index=False)
    CON.execute("CREATE INDEX med_case ON medication(CaseID,Offset)")
    CON.execute("CREATE INDEX range_case ON data_range(CaseID,Offset)")
    CON.execute("CREATE INDEX unit_case ON unitlog(CaseID,Offset)")
    CON.commit()
    vaso=m[m.DrugID.isin({1502,1550,1562,1593,1618})].copy()
    vaso=vaso.merge(valid[["CaseID","ICUOffset","TimeOfStay","WeightOnAdmission"]],on="CaseID")
    q={"source_medication_rows":medrows,"retained_medication_rows":len(m),
       "eligible_case_boundaries_age_weight":len(valid),"five_drug_rows":len(vaso),
       "invalid_interval_rows":int((vaso.OffsetDrugEnd<=vaso.Offset).sum()),
       "missing_interval_rows":int((vaso.OffsetDrugEnd.isna()|vaso.Offset.isna()).sum()),
       "nonpositive_rate_rows":int((vaso.AmountPerMinute<=0).sum()),
       "missing_rate_rows":int(vaso.AmountPerMinute.isna().sum())}
    v=vaso[(vaso.OffsetDrugEnd>vaso.Offset)&(vaso.AmountPerMinute>0)&(vaso.Offset<vaso.TimeOfStay)&(vaso.OffsetDrugEnd>vaso.ICUOffset)].copy()
    v["start"]=v[["Offset","ICUOffset"]].max(axis=1)
    v["end"]=v[["OffsetDrugEnd","TimeOfStay"]].min(axis=1)
    v=v.sort_values(["CaseID","DrugID","start","end"])
    prev=v.groupby(["CaseID","DrugID"]).end.transform(lambda x:x.cummax().shift())
    overlaps=v.start<prev
    q["same_drug_overlap_rows"]=int(overlaps.sum())
    q["same_drug_overlap_cases"]=int(v.loc[overlaps,"CaseID"].nunique())
    q["positive_icu_case_intersecting_drug_rows"]=len(v)
    q["single_dose_rows_in_icu_case"]=int((v.IsSingleDose==1).sum())
    q["single_dose_duration_60s_rows"]=int(((v.IsSingleDose==1)&((v.OffsetDrugEnd-v.Offset)==60)).sum())
    v.loc[overlaps,["CaseID","DrugID","start","end"]].to_sql("same_drug_overlaps",CON,index=False)
    source_ids=set(v.CaseID.astype(int))
    pd.DataFrame({"CaseID":sorted(source_ids)}).to_sql("signal_target_cases",CON,index=False)
    q["signal_target_cases"]=len(source_ids)
    CON.commit()
    export_json("medication_basic_checks.json",q)
    log(json.dumps(q))

def stream_signals():
    if already("signals"):
        raise RuntimeError("Signal output exists; do not append a second scan")
    ids=set(pd.read_sql_query("SELECT CaseID FROM signal_target_cases",CON).CaseID.astype(int))
    dictionary=pd.read_sql_query("SELECT * FROM dictionary",CON)
    counts=Counter()
    nonnull=Counter()
    retained=Counter()
    total=0
    selected=0
    # All IDs are counted, even if not selected for a candidate's feature window.
    # rawdata is deliberately not selected, decoded or persisted.
    cols=["CaseID","DataID","Offset","Val","cnt"]
    for i,c in enumerate(read("data_float_h",usecols=cols,chunksize=CHUNK)):
        total+=len(c)
        counts.update(c.DataID.value_counts().to_dict())
        nonnull.update(c.loc[c.Val.notna(),"DataID"].value_counts().to_dict())
        k=c[c.CaseID.isin(ids)&c.DataID.isin(SIGNALS)]
        if len(k):
            k.to_sql("signals",CON,index=False,if_exists="append",chunksize=10000)
            selected+=len(k)
            retained.update(k.DataID.value_counts().to_dict())
        if i%20==0:
            export_json("signal_scan_progress.json",{"status":"RUNNING","rows_scanned":total,"rows_selected":selected})
            log(f"SIGNALS rows={total} selected={selected}")
    CON.execute("CREATE INDEX signal_case_time ON signals(CaseID,Offset,DataID)")
    CON.commit()
    out=pd.DataFrame([{"DataID":int(k),"all_rows":int(n),"nonnull_Val":int(nonnull[k]),"selected_rows":int(retained[k])} for k,n in sorted(counts.items())])
    out=out.merge(dictionary[["ReferenceGlobalID","ReferenceValue","ReferenceName","ReferenceUnit"]],left_on="DataID",right_on="ReferenceGlobalID",how="left")
    out.to_csv(OUT/"signal_item_coverage.csv",index=False)
    export_json("signal_scan_progress.json",{"status":"COMPLETE","rows_scanned":total,"rows_selected":selected,"rawdata_decoded":False})
    log(f"SIGNALS_COMPLETE rows={total} selected={selected}")

def stream_labs():
    if already("labs"):
        raise RuntimeError("Lab output already exists")
    ids=set(pd.read_sql_query("SELECT CaseID FROM signal_target_cases",CON).CaseID.astype(int))
    counts=Counter()
    total=selected=0
    for i,c in enumerate(read("laboratory",chunksize=CHUNK)):
        total+=len(c)
        counts.update(c.LaboratoryID.value_counts().to_dict())
        k=c[c.CaseID.isin(ids)&c.LaboratoryID.isin(LABS)]
        if len(k):
            k.to_sql("labs",CON,index=False,if_exists="append",chunksize=10000)
            selected+=len(k)
    CON.execute("CREATE INDEX lab_case_time ON labs(CaseID,Offset,LaboratoryID)")
    CON.commit()
    pd.DataFrame([{"LaboratoryID":int(k),"rows":int(n)} for k,n in sorted(counts.items())]).to_csv(OUT/"laboratory_item_coverage.csv",index=False)
    export_json("laboratory_scan_summary.json",{"status":"COMPLETE","rows_scanned":total,"rows_selected":selected})
    log(f"LABS_COMPLETE rows={total} selected={selected}")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=["prep","signals","labs"])
    a=p.parse_args()
    {"prep":prep,"signals":stream_signals,"labs":stream_labs}[a.phase]()
    CON.close()
