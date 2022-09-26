import os

# Only need to import the modules if NOT running under PyTest
# from: https://stackoverflow.com/a/58866220/13885200
if "PYTEST_CURRENT_TEST" not in os.environ:
    import database_convert
    import input_database
    import output_database
    import taxon_ids
