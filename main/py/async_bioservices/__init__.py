import sys

# Only need to import the modules if NOT running under PyTest
# from: https://stackoverflow.com/questions/25188119
if "pytest" not in sys.modules:
    import database_convert
    import input_database
    import output_database
    import taxon_ids
else:
    from . import database_convert
    from . import input_database
    from . import output_database
    from . import taxon_ids
