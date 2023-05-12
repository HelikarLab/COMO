import io
import sys
import logging
import contextlib
from typing import Iterator


def stringlist_to_list(stringlist: str | list[str]) -> list[str]:
    """
    We are attempting to move to a new method of gathering a list of items from the command line
    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list, assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    if isinstance(stringlist, str):
        new_list: list[str]
        if stringlist.startswith("[") and stringlist.endswith("]"):
            # Remove any brackets from the first and last items; replace quotation marks and commas with nothing
            new_list = stringlist.strip("[]").replace("'", "").replace(" ", "").split(",")
            
            # Show a warning if more than one item is present in the list (this means we are using the old method)
            print("DeprecationWarning: Please use the new method of providing context names, i.e. --output-filetypes 'type1 type2 type3'.")
            print("If you are using COMO, this can be done by setting the 'context_names' variable to a simple string separated by spaces. Here are a few examples!")
            print("context_names = 'cellType1 cellType2 cellType3'")
            print("output_filetypes = 'output1 output2 output3'")
            print("\nYour current method of passing context names will be removed in the future. Please update your variables above accordingly!\n\n")
            
        else:
            new_list = stringlist.split(" ")

        return new_list
    return stringlist

@contextlib.contextmanager
def suppress_stdout() -> Iterator[None]:
    with io.StringIO() as buffer:
        try:
            sys.stdout = buffer
            yield
        
        finally:
            sys.stdout = sys.__stdout__
