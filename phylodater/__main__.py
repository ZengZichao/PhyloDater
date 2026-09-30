"""
PhyloDater - 多软件并行系统发育定年平台

允许用户通过 python -m phylodater 运行
"""

import sys

from phylodater.cli.main import main

if __name__ == "__main__":
    sys.exit(main())
