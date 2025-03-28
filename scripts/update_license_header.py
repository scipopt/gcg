#!/usr/bin/env python3

import os
import re
import sys
import logging
import argparse


# #############################################################################
#                                   Constants
# #############################################################################

DEFAULT_LOG_LVL         = logging.INFO
DEFAULT_LOGFILE_NAME    = "update_license_header.log"
DEFAULT_LOG_TO_STDOUT   = True

DEFAULT_PATHS           = ["check", "doc", "src", "tests", "stats"]
DEFAULT_FILE_EXTENSIONS = ["c", "h", "cpp", "hpp", "awk", "sh", "dag", ""]

DEFAULT_LOG_PADDING     = 28
DEFAULT_DRY_RUN         = False
DEFAULT_IGNORE_ACT_LIC  = False

REGEX_LIC_HEADER_C      = re.compile(r"^\/(\*\s){38}\*.*\n(^\/\*\s[^\*].*\n)*^\/(\*\s){38}\*.*", 
                                     flags=re.MULTILINE)
REGEX_LIC_HEADER_BASH   = re.compile(r"^\#(\*\s){38}\*.*\n(^\#\*\s[^\*].*\n)*^\#(\*\s){38}\*.*", 
                                     flags=re.MULTILINE)


# #############################################################################
#                                   Functions
# #############################################################################

def setup_logging(logfile, log_to_stdout, log_level):
    """ Sets up logging to a file and optionally the console. """
    log_format = "%(asctime)s - [ %(levelname)8s ] - %(message)s "
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        format=log_format,
        level=log_level
    )

    if log_to_stdout:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(log_format))
        logging.getLogger().addHandler(console_handler)


def span_to_lines(string, span):
    """ Converts a span in a string to a span in lines. """
    line_start = string.count('\n', 0, span[0]) + 1
    line_end = string.count('\n', 0, span[1]) + 1
    return (line_start, line_end)


def replace_license_header(path_to_gcg_root, paths_from_gcg_root, file_ext_to_consider, path_to_new_lic_header, dry_run, ignore_act_lic):
    """ Replaces the license header in the files in the given paths. """
    logging.info(f"*** Replacing license in paths: {paths_from_gcg_root} ***")
    logging.info(f"*** Replacing license in file extensions: {file_ext_to_consider} ***")

    match_no_ext = ("" in file_ext_to_consider)
    if match_no_ext:
        file_ext_to_consider.remove("")
    
    # Read new license header
    with open(path_to_new_lic_header, "r") as f:
        new_license_header_c = f.read()
    new_license_header_c = new_license_header_c.rstrip('\n')
    new_license_header_bash = re.sub(r"^\/", r"#", new_license_header_c)
    new_license_header_bash = re.sub(r"\/\n\/", r"\n#", new_license_header_bash, 
                                     flags=re.MULTILINE)
    new_license_header_bash = re.sub(r"\/$", r"", new_license_header_bash)

    match_replace_pairs = [
        (REGEX_LIC_HEADER_C, new_license_header_c),
        (REGEX_LIC_HEADER_BASH, new_license_header_bash)
        ]

    # Find and replace license header
    for path in paths_from_gcg_root:
        logging.info(f"*** Processing path: {path} ***")
        for root, dirs, files in os.walk(os.path.join(path_to_gcg_root, path)):
            for file in files:
                name, ext = os.path.splitext(file)
                ext = ext.strip(".")
                if ext in file_ext_to_consider or (match_no_ext and not "." in file):
                    file_path = os.path.join(root, file)
                    logging.info("Checking file: ".ljust(DEFAULT_LOG_PADDING) + f"{file_path}")
                    try:
                        with open(file_path, "r", encoding="utf-8") as f:
                            content = f.read()
                    except UnicodeDecodeError:
                        # If file is not UTF-8, skip it
                        logging.warning("Skipping non-UTF-8 file: ".ljust(DEFAULT_LOG_PADDING) + f"{file_path}")
                        continue
                    # Match regex for license header
                    for regex, new_lic_header in match_replace_pairs:
                        result = re.search(regex, content)
                        if result:
                            if ignore_act_lic and result.group() == new_lic_header:
                                logging.info("License up to date. Ignore: ".ljust(DEFAULT_LOG_PADDING) + f"{file_path}")
                                continue
                            logging.info("Updating license in: ".ljust(DEFAULT_LOG_PADDING) + f"{file_path}")
                            l1, l2 = span_to_lines(content, result.span())
                            logging.debug("Match lines: ".ljust(DEFAULT_LOG_PADDING) + f"{l1}-{l2}")
                            logging.debug(f"Replacing matched header:\n{result.group()}")
                            logging.debug(f"Replacing with:\n{new_lic_header}")
                            if not dry_run:
                                content = re.sub(regex, new_lic_header, content)
                                with open(os.path.join(root, file), "w") as f:
                                    f.write(content)
                            break


# #############################################################################
#                                   Main
# #############################################################################

def parse_arguments():
    """ Parses command-line arguments. """
    parser = argparse.ArgumentParser(description="Replace license headers in GCG project files.")

    parser.add_argument("path_to_gcg_root", type=str, 
                        help="Path to the GCG project root.")
    parser.add_argument("path_to_new_license", type=str, 
                        help="Path to the new license header file. (I.e., C-comment style formatted license header.)")

    parser.add_argument("--logfile", type=str, default=DEFAULT_LOGFILE_NAME, 
                        help="Path to the log file.")
    parser.add_argument("--loglevel", type=str, choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                        default="INFO", 
                        help="Logging level. Default: INFO. (Set to CRITICAL to disable logging.)")
    parser.add_argument("--log-to-stdout", default=DEFAULT_LOG_TO_STDOUT, action="store_true", 
                        help="Enable logging to stdout.")

    parser.add_argument("--paths", nargs="+", default=DEFAULT_PATHS, 
                        help=f"Directories (paths relative to GCG root) to search for files. (Default: {DEFAULT_PATHS}).")
    parser.add_argument("--extensions", nargs="+", default=DEFAULT_FILE_EXTENSIONS, 
                        help=f"File extensions to consider (e.g., c h cpp). Use '' to match " +\
                              "files with no extension. (Default: {DEFAULT_FILE_EXTENSIONS}).")
    parser.add_argument("--dry-run", default=DEFAULT_DRY_RUN, action="store_true", 
                        help="Do not update license headers. Only find them.")
    parser.add_argument("--ignore-act-lic", default=DEFAULT_IGNORE_ACT_LIC, action="store_true", 
                        help="Ignore license headers that match the one in the act-lic file.")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    # Validate paths
    if not os.path.isdir(args.path_to_gcg_root):
        sys.exit(f"Error: GCG root directory not found: {args.path_to_gcg_root}")
    if not os.path.isfile(args.path_to_new_license):
        sys.exit(f"Error: License header file not found: {args.path_to_new_license}")

    # Correct log level in case find-not-update is set
    # (In that case, we only want to find license headers, but not update them.)
    # (So we also want to see differences in the license headers logged as warnings.)
    if args.dry_run and args.loglevel == "INFO":
        args.loglevel = "DEBUG"

    # Setup logging
    setup_logging(args.logfile, args.log_to_stdout, getattr(logging, args.loglevel))

    # Run license replacement
    replace_license_header(args.path_to_gcg_root, 
                           args.paths, 
                           args.extensions, 
                           args.path_to_new_license, 
                           args.dry_run,
                           args.ignore_act_lic)

