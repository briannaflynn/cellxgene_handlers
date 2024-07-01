import logging
import sys
import io
from contextlib import contextmanager
from datetime import datetime

logger_ = logging.getLogger("cxghandler")
# check if logger has been initialized
if not logger_.hasHandlers() or len(logger_.handlers) == 0:
    logger_.propagate = False
    logger_.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(name)s - %(levelname)s - %(message)s", datefmt="%H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger_.addHandler(handler)

def setup_logger(name, log_file, level=logging.DEBUG):
    """Function to setup as many loggers as you want"""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers = []  # Clear existing handlers

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

class StreamToLogger(object):
    """
    Context manager for redirecting stdout to a logger.
    """
    def __init__(self, logger, level=logging.INFO):
        self.logger = logger
        self.level = level
        self.linebuf = ''

    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self  # Redirect stdout to this class instance

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.stdout  # Reset stdout to its original value

    def write(self, message):
        if message.rstrip() != "":
            self.logger.log(self.level, message.rstrip())

    def flush(self):
        pass


