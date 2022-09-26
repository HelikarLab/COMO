import sys

# Only need to import the modules if NOT running under PyTest
# from: https://stackoverflow.com/questions/25188119
if "pytest" not in sys.modules:
    import Crux
    import FileInformation
    import FTPManager
    import proteomics_preprocess
