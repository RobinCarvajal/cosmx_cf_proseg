
import sqlite3
from pathlib import Path

# defien main dir
DATA_DIR = Path('/Volumes/robin_work/cosmx_gray/data')
DB_PATH = DATA_DIR / 'cf.db'
# define connections and cursor
conn = sqlite3.connect(DB_PATH)
cur = conn.cursor()

# your dict
epithelial_markers = {
    'undefined': ['COX1', 'COX2'],
    'epithelial': ['EPCAM', 'KRT14', 'CDH1', 'CD24'],
    'basal': ['KRT5', 'KRT14', 'KRT15', 'KRT17', 'IL33'],
    'secretory': ['SCGB1A1', 'SCGB3A1'],
    'ciliated': ['FOXJ1', 'TUBB4B', 'SAA1'],
    'AT1': ['RTKN2', 'COL4A3', 'FSTL3', 'AGER'], 
    'AT2': ['NKX2-1', 'ETV5', 'NAPSA', 'LAMP3', 'CD36', 'LPCAT1', 'SFTPD'],
}

# create table if it doesn't exist
cur.execute("""
CREATE TABLE IF NOT EXISTS epithelial (
    celltype TEXT PRIMARY KEY,
    markers TEXT
)
""")

# insert data
for celltype, markers in epithelial_markers.items():
    markers_str = "/".join(markers)
    cur.execute(
        "INSERT OR REPLACE INTO epithelial (celltype, markers) VALUES (?, ?)",
        (celltype, markers_str)
    )

# save changes and close
conn.commit()
conn.close()
