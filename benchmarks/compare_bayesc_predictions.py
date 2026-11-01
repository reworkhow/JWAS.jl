"""Compare the same BayesC fixtures from upstream and candidate projects."""
import csv,json,math,sys
from pathlib import Path
left,right,report=map(Path,sys.argv[1:])
report.mkdir(parents=True,exist_ok=True)
def keyed(path,keys,expected):
    rows=list(csv.DictReader(path.open()))
    data={}
    for r in rows:
        key=tuple(r[k].strip().strip('"') for k in keys)
        assert key not in data,(path,"duplicate key",key)
        data[key]=r
    assert len(rows)==len(expected) and set(data)==expected,(path,"wrong keys/count")
    return data
cases=list(csv.DictReader((left/"cases.csv").open()))
assert cases==list(csv.DictReader((right/"cases.csv").open()))
results=[]
for case in cases:
    name=case["name"];traits=int(case["traits"])
    result=dict(case)
    for kind,keys,expected,columns in (
        ("animals",("Trait","ID"),{(f"y{k}",f"animal_{i:04d}") for k in range(1,traits+1) for i in range(29,41)},("EBV","PEV")),
        ("markers",("Trait","Marker_ID"),{(f"y{k}",f"rs_{1000+7*i}") for k in range(1,traits+1) for i in range(1,20)},("Estimate","SD","Model_Frequency")),
    ):
        a=keyed(left/f"{name}_{kind}.csv",keys,expected)
        b=keyed(right/f"{name}_{kind}.csv",keys,expected)
        assert a.keys()==b.keys()
        aligned=[]
        for key in sorted(expected):
            row=dict(zip(keys,key))
            for col in columns:
                u,v=float(a[key][col]),float(b[key][col])
                assert math.isfinite(u) and math.isfinite(v)
                row[f"upstream_{col}"]=u;row[f"candidate_{col}"]=v
            aligned.append(row)
        for col in columns:
            delta=[abs(r[f"upstream_{col}"]-r[f"candidate_{col}"]) for r in aligned]
            result[f"max_abs_{col}"]=max(delta)
            result[f"exact_{col}"]=all(d==0 for d in delta)
            result[f"match_{col}"]=all(math.isclose(r[f"upstream_{col}"],r[f"candidate_{col}"],abs_tol=1e-4,rel_tol=1e-4) for r in aligned)
        with (report/f"{name}_{kind}.csv").open("w") as handle:
            w=csv.DictWriter(handle,fieldnames=list(aligned[0]),lineterminator="\n")
            w.writeheader();w.writerows(aligned)
        if case is cases[0]:print(kind,"joined sample",aligned[:2])
    results.append(result)
with (report/"summary.csv").open("w") as handle:
    w=csv.DictWriter(handle,fieldnames=list(results[0]),lineterminator="\n")
    w.writeheader();w.writerows(results)
cols=("EBV","PEV","Estimate","SD","Model_Frequency")
summary={"cases":len(results),"all_pass":all(r[f"match_{c}"] for r in results for c in cols),
    "all_exact":all(r[f"exact_{c}"] for r in results for c in cols),
    "max_abs":{c:max(r[f"max_abs_{c}"] for r in results) for c in cols},
    "failures":[r["name"] for r in results if not all(r[f"match_{c}"] for c in cols)]}
(report/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
assert summary["all_pass"]
