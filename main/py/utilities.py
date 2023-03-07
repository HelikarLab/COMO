def stringlist_to_list(stringlist: list[str]) -> list[str]:
    """
    We are attempting to move to a new method of gathering a list of items from the command line
    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list, assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    first_context: str = stringlist[0]
    last_context: str = stringlist[-1]
    
    if first_context.startswith("[") and last_context.endswith("]"):
        print("DeprecationWarning: Please use the new method of providing context names, i.e. --output-filetypes 'typ1 type 2 type3'.")
        print("If you are using COMO, this can be done by setting the 'context_names' variable to a simple string separated by spaces. Here are a few examples!")
        print("context_names = 'cellType1 cellType2 cellType3'")
        print("output_filetypes = 'output1 output2 output3'")
        print("\nYour current method of passing context names will be removed in the future. Please update your variables above accordingly!\n\n")
        
        # Remove any brackets from the first and last items; replace quotation marks and commas with nothing
        return [
            i.strip("[]").replace("'", "").replace(",", "")
            for i in stringlist
        ]
    return stringlist
