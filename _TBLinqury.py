import sqlite3

def print_table_headers(db_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get all table names
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()

    if not tables:
        print("No tables found in the database.")
        return

    print(f"\nDatabase: {db_path}")
    print("=" * (len(db_path) + 10))

    # Loop through each table and print its columns
    for table_name_tuple in tables:
        table_name = table_name_tuple[0]
        print(f"\nTable: {table_name}")
        print("-" * (len(table_name) + 8))

        # Get the column info
        cursor.execute(f"PRAGMA table_info({table_name});")
        columns = cursor.fetchall()

        if not columns:
            print("  (No columns found)")
        else:
            # Print column headers
            print("  Columns:")
            for col in columns:
                # col = (cid, name, type, notnull, dflt_value, pk)
                print(f"    - {col[1]} ({col[2]})")

    # Clean up
    conn.close()

if __name__ == "__main__":
    # db_file = "Ligases/eliah.db"
    db_file = "Ligases/Ligase_Recruiter.db"
    # db_file = "Ligases/Ligase_Recruiter_2.db"

    print_table_headers(db_file)







