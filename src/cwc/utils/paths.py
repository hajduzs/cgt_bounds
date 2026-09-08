import os
from pathlib import Path

# Base directory of the project (the repository root)
BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data"

CODES_DIR = DATA_DIR / "codes"
BLOCKCODES_DIR = DATA_DIR / "blockcodes"
CODE_LENGTHS_JSON = DATA_DIR / "code_lengths.json"
BLOCKCODE_LENGTHS_JSON = DATA_DIR / "blockcode_lengths.json"
LINEAR_QUARY_CODES_JSON = DATA_DIR / "linear_quary_codes.json"
TAPESTRY_CODES_JSON = DATA_DIR / "codes_from_tapestry2.json"
PARTITIONS_DIR = DATA_DIR / "partitions"

def ensure_dirs():
    os.makedirs(CODES_DIR, exist_ok=True)
    os.makedirs(CODES_DIR / "old", exist_ok=True)
    os.makedirs(CODES_DIR / "subopt", exist_ok=True)
    os.makedirs(BLOCKCODES_DIR, exist_ok=True)
    os.makedirs(BLOCKCODES_DIR / "subopt", exist_ok=True)
    os.makedirs(PARTITIONS_DIR, exist_ok=True)
