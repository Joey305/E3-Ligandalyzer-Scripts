#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build_Ligase_Database.py
----------------------------------------
Creates a unified SQLite database (Ligase_Recruiter.db)
from all CSV files in Ligase_Table/.

Each CSV becomes its own table named after the file.
"""

import sqlite3
import pandas as pd
from pathlib import Path

# ============================================================
# ⚙️ Configuration
# ============================================================
TABLE_DIR = Path("Ligase_Table")
DB_PATH = Path("Ligase_Recruiter.db")
REPLACE_EXISTING = True  # overwrite DB if it already exists

# ============================================================
# 🚀 Build Database
# ============================================================
def main():
    if not TABLE_DIR.exists():
        print(f"❌ Directory not found: {TABLE_DIR.resolve()}")
        return

    if DB_PATH.exists() and REPLACE_EXISTING:
        DB_PATH.unlink()
        print(f"🗑️ Removed existing {DB_PATH.name}")

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    csv_files = sorted(TABLE_DIR.glob("*.csv"))
    if not csv_files:
        print("⚠️ No CSV files found in Ligase_Table/")
        return

    print(f"📂 Found {len(csv_files)} CSV files to load into {DB_PATH.name}\n{'='*70}")

    for f in csv_files:
        table_name = f.stem
        try:
            df = pd.read_csv(f)
            df.to_sql(table_name, conn, if_exists="replace", index=False)
            print(f"✅ Loaded {table_name:<35} → {len(df):>6,} rows × {len(df.columns):>3} cols")
        except Exception as e:
            print(f"❌ Failed to load {table_name}: {e}")

    # Print final summary of tables
    print("\n📜 Database Summary:")
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    for (tbl,) in cursor.fetchall():
        cursor.execute(f"SELECT COUNT(*) FROM {tbl}")
        (count,) = cursor.fetchone()
        print(f"   • {tbl:<35} {count:>6,} rows")

    conn.close()
    print(f"\n✅ Database built successfully → {DB_PATH.resolve()}")

if __name__ == "__main__":
    main()
