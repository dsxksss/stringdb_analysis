import logging
from .utils import set_logging_default_config


def main() -> int:
    set_logging_default_config()
    logging.info("My Python Script Template!")
    return 0
