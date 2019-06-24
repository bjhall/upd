#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
upd.__main__
~~~~~~~~~~~~~~~~~~~~~

The main entry point for the command line interface.

Invoke as `upd` (if installed)
or ``python -m upd`` (no install required).
"""
import sys

from upd.cli import cli as base_command


if __name__ == '__main__':
    # exit using whatever exit code the CLI returned
    sys.exit(base_command())
