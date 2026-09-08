"""Put src/ on sys.path so tests import the pipeline the same way notebooks do.

`src/` uses flat imports (`from qubits import ...`) and is deliberately not a package,
so it must be on sys.path -- see CLAUDE.md.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
