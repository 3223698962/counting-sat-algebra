"""Allow `python -m counting_sat_algebra`."""

import sys

from .verify import main

if __name__ == "__main__":
    sys.exit(main())
