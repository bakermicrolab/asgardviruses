#!/usr/bin/python3

"""
Motivation: functions for command-line arguments.

Author: Ian Rambo

Contact: ian.rambo@utexas.edu, imrambo@lbl.gov
"""

import magic
import gzip
import os
import logging
import re

def optstring_join(optdict):
    """
    Join a dictionary of command line options into a single string.
    """
    optstring = ' '.join([str(param) + ' ' + str(val) for param, val in optdict.items()])
    return optstring
#------------------------------------------------------------------------------
def exec_cmd_generate(exec_path, optdict):
    """
    Create an argument list from a dictonary to use with subprocess.run()
    For use when running an executable with shell=False
    """
    optlist = []
    for param, val in optdict.items():
        optlist.append(str(param))
        optlist.append(str(val))
    optlist.insert(0, exec_path)
    return optlist
#------------------------------------------------------------------------------
def is_gzipped(file_path):
    """
    Test if a file is gzipped.
    """
    is_gzip = False
    if magic.from_file(file_path).startswith('gzip compressed data'):
        is_gzip = True
        return is_gzip
    else:
        return is_gzip
#------------------------------------------------------------------------------
def joblog_test(joblog):
    """
    Parse GNU Parallel joblog to recover commands with non-zero exit codes.
    Returns a two list of tuples containing the command and exit code
    (zero and non-zero exit codes).
    """
    nonzero = []
    zero = []
    if os.path.exists(joblog) and os.path.getsize(joblog) > 0:
        with open(joblog, 'r') as jl:
            next(jl)
            for record in jl:
                record_list = record.split()
                command = record_list[-1]
                exit_code = record_list[6]
                if exitCode != 0:
                    nonzero.append((command, exit_code))
                else:
                    zero.append([command, exit_code])
    # if nonzero:
    #     warning_message = '%d of %d commands returned a non-zero exit code' % (len(nonzero), len(nonzero)+len(zero))
    #     logging.warning(warning_message)
    #     for c,e in nonzero:
    #         warning_message = 'Command: "%s" returned a non-zero exit code of: %s' % (c, e)
    #         logging.warning(warning_message)
    return nonzero, zero
#------------------------------------------------------------------------------
def get_basename(file_path):
    try:
        basename = os.path.basename(file_path)
        #Remove two extensions, e.g. foo.tar.gz becomes foo
        if re.match(r'^.*?\.[a-z]+\.[a-z]+$', basename):
            basename = re.findall(r'^(.*?)\.[a-z]+\.[a-z]+$', basename)[0]
        else:
            basename = os.path.splitext(basename)[0]
    except:
        print('no file path specified!')
    return basename
