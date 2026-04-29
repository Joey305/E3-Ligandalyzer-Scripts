import pandas as pd
from pathlib import Path

inp = Path("Ligase_Table/Ligase_Recruiters_Superclustered.csv")
out = Path("Ligase_Table/Scaffold_Unified_Map.csv")

df = pd.read_csv(inp)

required = {"Unified_Scaffold_ID", "Scaffold_Hash", "Scaffold_ID"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns from {inp}: {missing}")

unified = (
    df[["Unified_Scaffold_ID", "Scaffold_Hash", "Scaffold_ID"]]
    .drop_duplicates()
    .rename(columns={"Scaffold_ID": "Original_Scaffold_ID"})
    .sort_values(["Unified_Scaffold_ID", "Original_Scaffold_ID"])
    .reset_index(drop=True)
)

unified.to_csv(out, index=False)

print(f"✅ Wrote {len(unified)} rows to {out}")
print(f"✅ Unique unified scaffold IDs: {unified['Unified_Scaffold_ID'].nunique()}")
print(f"✅ Unique scaffold hashes: {unified['Scaffold_Hash'].nunique()}")

bad = (
    unified.groupby("Scaffold_Hash")["Unified_Scaffold_ID"]
    .nunique()
    .reset_index(name="Num_Unified_IDs")
)
bad = bad[bad["Num_Unified_IDs"] > 1]

if len(bad):
    print("❌ Some hashes map to multiple unified scaffold IDs:")
    print(bad.head(20).to_string(index=False))
else:
    print("🎯 PASS: each Scaffold_Hash maps to exactly one Unified_Scaffold_ID")
