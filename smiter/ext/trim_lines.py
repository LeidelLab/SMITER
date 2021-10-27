#!/usr/bin/env python
import errno
import fileinput
import sys


def main():
    try:
        for i, line in enumerate(fileinput.input()):
            if i == 0:
                sys.stdout.write(line)
                continue
            line = line.split(",")[2:]
            sys.stdout.write(",".join(line))
        sys.stdout.flush()
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(1)


if __name__ == "__main__":
    main()
