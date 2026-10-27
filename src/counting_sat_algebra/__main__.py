"""Allow `python -m counting_sat_algebra`."""

import sys

def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] == 'count':
        from .cli import count_main
        return count_main(args[1:])
    from .verify import main as verify_main
    return verify_main(args[1:] if args and args[0] == 'verify' else args)

if __name__ == "__main__":
    sys.exit(main())
