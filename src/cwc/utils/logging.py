import logging
import sys

class CWCLogger(logging.Logger):
    def _log_with_args(self, level, args, kwargs, indent=""):
        if args:
            msg = indent + " ".join(map(str, args))
            super()._log(level, msg, (), **kwargs)

    def info(self, *args, **kwargs):
        if self.isEnabledFor(logging.INFO):
            self._log_with_args(logging.INFO, args, kwargs)

    def debug(self, *args, **kwargs):
        if self.isEnabledFor(logging.DEBUG):
            self._log_with_args(logging.DEBUG, args, kwargs)

    def warning(self, *args, **kwargs):
        if self.isEnabledFor(logging.WARNING):
            self._log_with_args(logging.WARNING, args, kwargs)

    def error(self, *args, **kwargs):
        if self.isEnabledFor(logging.ERROR):
            self._log_with_args(logging.ERROR, args, kwargs)

# Register the custom logger class
logging.setLoggerClass(CWCLogger)

# Initialize the 'cwc' logger
logger = logging.getLogger("cwc")

# Configure default logging to stdout if no handlers are present
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # Default level matches original debug=4 (log1, log2, log3 visible)
    logger.setLevel(logging.INFO)

# Compatibility functions
def log1(*args):
    logger.info(*args)

def log2(*args):
    # We add the indent here to match original behavior
    logger.info(" ", *args)

def log3(*args):
    logger.info("  ", *args)

def log4(*args):
    logger.debug("   ", *args)

def log5(*args):
    logger.debug("    ", *args)

def log6(*args):
    logger.debug("     ", *args)

def is_debug(level):
    if level <= 3:
        return logger.getEffectiveLevel() <= logging.INFO
    else:
        return logger.getEffectiveLevel() <= logging.DEBUG

def set_debug(level):
    if level > 3:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

def inc_debug():
    logger.setLevel(logging.DEBUG)

def dec_debug():
    logger.setLevel(logging.INFO)
