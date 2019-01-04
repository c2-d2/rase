#! /usr/bin/env python3

"""rase

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import re
import sys

sys.path.append(os.path.dirname(__file__))
import version

PROGRAM='rase'
VERSION=version.VERSION
DESC=''

def parse_args():

    class CustomArgumentParser (argparse.ArgumentParser):
        def print_help(self):
            msg=self.format_help()
            msg=msg.replace("usage:", "Usage:  ")
            for x in 'PY_EXPR', 'PY_CODE':
                msg=msg.replace("[{x} [{x} ...]]\n            ".format(x=x), x)
                msg=msg.replace("[{x} [{x} ...]]".format(x=x), x)
            repl=re.compile(r'\]\s+\[')
            msg=repl.sub("] [",msg)
            msg=msg.replace("\n  -0","\n\nAdvanced options:\n  -0")
            msg=msg.replace(" [-h] [-v]","")
            msg=msg.replace("[-0","\n                    [-0")
            print(msg)

        def format_help(self):
            formatter = self._get_formatter()
            formatter.add_text(" \n"+self.description)
            formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
            formatter.add_text(self.epilog)

            # positionals, optionals and user-defined groups
            for action_group in self._action_groups:
                formatter.start_section(action_group.title)
                formatter.add_text(action_group.description)
                formatter.add_arguments(action_group._group_actions)
                formatter.end_section()

            return formatter.format_help()

    parser = CustomArgumentParser (
            formatter_class=argparse.RawTextHelpFormatter,
            description=
            "Program: {} ({})\n".format(PROGRAM, DESC)+
            "Version: {}\n".format(VERSION) +
            "Author:  Karel Brinda <kbrinda@hsph.harvard.edu>",
            )
    parser.add_argument('-v', '--version',
            action='version',
            version='{} {}'.format(PROGRAM, VERSION),
            )

    args = parser.parse_args()

    return args


def main ():
    args=parse_args()
    print("Hello")


if __name__ == "__main__":
    main()
