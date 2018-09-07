import os
import platform
import json
import argparse
import logging

from .worker.worker import GlycomodWorker
from .worker.utils import check_internet_conn
from .worker.utils import validate_filename
from .worker.utils import check_text_input_consistency 

logging.basicConfig(format='[%(levelname)s][%(funcName)s] %(message)s', level=logging.INFO)
logger = logging.getLogger(__file__)

def init_config(path=None):
    if not path:
        conf_path = os.path.join(os.path.dirname(__file__), "worker","config.json")
    elif path.endswith(".json"):
        conf_path = os.path.abspath(path)
    else:
        raise EnvironmentError(f"Initializing configuration failed with file: {path}")
    with open(conf_path, 'r') as config:
        cfg = json.load(config)
    logger.debug("Initialized config: %s:" % conf_path)
    return cfg

def which_driver():
    cur = os.path.dirname(os.path.abspath(__file__))
    system = platform.system()
    if system == "Linux":
        ch_drvr = os.path.normpath(os.path.join(cur, "worker", "dep", "linux", "chromedriver"))
    elif system == "Darwin":
        ch_drvr = os.path.normpath(os.path.join(cur, "worker", "dep", "mac", "chromedriver"))
    elif system == "Windows":
        ch_drvr = os.path.normpath(os.path.join(cur, "worker", "dep", "win", "chromedriver.exe"))
    else:
        raise EnvironmentError(f"Unsupported operating system: {system}")
    return ch_drvr

# TODO test main()
# TODO make GUI
# TODO setup.py
# TODO write readme.md
# TODO choose licese
def main():
    driver_path = which_driver()
    if check_internet_conn():
        parser = argparse.ArgumentParser(description='Provide a filepath and run Glycomod search and save results as text or table.')
        parser.add_argument("--echo", "-e", action="store_true", help="Display results in terminal/console")
        parser.add_argument("--debug", "-d", action="store_true", help="Display debug logs in console")
        parser.add_argument("path", type=str, help="Path to text file containig glycan masses. Masses must be floating point numbers separated by newlines")
        parser.add_argument("--tag", "-t", type=str, help=f"N-glycan reducing end tag NAME.\nSupported: 2-AB, ProA")
        parser.add_argument("--filename", "-f", type=str, help="Filename used for saving .csv or .txt")
        parser.add_argument("--config", "-cfg", type=str, help="Path to config. If not provided deafault is used")
        parser.add_argument("--txt", action="store_true", help="Save results as text")
        args = parser.parse_args()
        if args.debug:
            logging.getLogger().setLevel(logging.DEBUG)
        if args.path.endswith(".txt"):
            path = os.path.normpath(args.path)
            with open(path, 'r') as infile:
                text_input = infile.read()
            check_text_input_consistency(text_input)
        else:
            raise TypeError("Provided file is not a textfile. Please provide a proper text(.txt) file.")
        if args.filename:
            validate_filename(args.filename)
            logger.debug("Using filename: %s" % args.filename)

        cfg = init_config(path=args.config)
        gw = GlycomodWorker(cfg, driver_path=driver_path, text_input=text_input, reducing_end=args.tag, save_txt=args.txt)
        logger.debug(f"Running GlycomodWorker\nParams: {args}")
        gw.run()
        if args.echo:
            gw.output_text(to_std_out=True)
        # check if config was modified - run some sort of checksum
            # check if stored checksum == calculated checksum
        logger.debug("SEARCH FINISHED SUCCESSFULLY")
    else:
        raise EnvironmentError("Unable to connect to Glycomod. No internet connection.")
    
if __name__ == "__main__":
    main()