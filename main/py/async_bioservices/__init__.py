import sys

# Only need to import the modules if NOT running under PyTest
# from: https://stackoverflow.com/questions/25188119
if "pytest" not in sys.modules:
    from async_bioservices import database_convert
    from async_bioservices import input_database
    from async_bioservices import output_database
    from async_bioservices import taxon_ids
