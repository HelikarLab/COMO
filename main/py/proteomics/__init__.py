import os

# Only need to import the modules if NOT running under PyTest
# from: https://stackoverflow.com/a/58866220/13885200
if "PYTEST_CURRENT_TEST" not in os.environ:
    import Crux
    import FileInformation
    import FTPManager
    import proteomics_preprocess
