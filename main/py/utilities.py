def stringlist_to_list(stringlist: str) -> list[str]:
    """
    We are attempting to move to a new method of gathering a list of items from the command line
    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list, assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    if "[" == stringlist[0] and "]" == stringlist[-1]:
        print("DEPRECATED: Please use the new method of providing context names, i.e. --output-filetypes 'typ1 type 2 type3'.")
        print("If you are using MADRID, this can be done by setting the 'context_names' variable to a simple string separated by spaces: output_filetypes='type1 type2 type3'")
        print("Your current method of passing context names will be removed in the future. Please update your variables above accordingly!\n\n")
        new_list = stringlist.rstrip("]").lstrip("[").replace("'", "").split(",")
        new_list = [x.strip() for x in new_list]
    else:
        new_list = stringlist.split(" ")

    return new_list
