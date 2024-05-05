from pathlib import Path
import sys


bin = Path(__file__).parent.parent.parent.joinpath("bin")
sys.path.append(str(bin))

import view_alignments
