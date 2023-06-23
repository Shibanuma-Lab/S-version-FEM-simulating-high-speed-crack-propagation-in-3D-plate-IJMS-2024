import logging
from logzero import setup_logger

logger = setup_logger(
    name='arrest-logger',
    logfile='logs/arrest.log',
    level=logging.INFO,
    formatter=None,
    fileLoglevel=logging.DEBUG,
    disableStderrLogger=False,
)
