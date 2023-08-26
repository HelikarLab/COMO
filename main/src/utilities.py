import io
import sys
import contextlib
import pandas as pd
from typing import Iterator


def stringlist_to_list(stringlist: str | list[str]) -> list[str]:
    """
    We are attempting to move to a new method of gathering a list of items from the command line
    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list, assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    if isinstance(stringlist, str):
        
        if stringlist.startswith("[") and stringlist.endswith("]"):
            # Remove any brackets from the first and last items; replace quotation marks and commas with nothing
            new_list: list[str] = stringlist.strip("[]").replace("'", "").replace(" ", "").split(",")
            
            # Show a warning if more than one item is present in the list (this means we are using the old method)
            print(
                "DeprecationWarning: Please use the new method of providing context names, i.e. --output-filetypes 'type1 type2 type3'.")
            print(
                "If you are using COMO, this can be done by setting the 'context_names' variable to a simple string separated by spaces. Here are a few examples!")
            print("context_names = 'cellType1 cellType2 cellType3'")
            print("output_filetypes = 'output1 output2 output3'")
            print(
                "\nYour current method of passing context names will be removed in the future. Please update your variables above accordingly!\n\n")
        
        else:
            new_list: list[str] = stringlist.split(" ")
        
        return new_list
    return stringlist


def split_gene_expression_data(expression_data, recon_algorithm="GIMME"):
    """
    Splits genes that have mapped to multiple Entrez IDs are formated as "gene12//gene2//gene3"
    """
    if recon_algorithm in ["IMAT", "TINIT"]:
        expression_data.rename(columns={"ENTREZ_GENE_ID": "Gene", "combine_z": "Data"}, inplace=True)
    else:
        expression_data.rename(columns={"ENTREZ_GENE_ID": "Gene", "Active": "Data"}, inplace=True)
    
    expression_data = expression_data.loc[:, ["Gene", "Data"]]
    expression_data["Gene"] = expression_data["Gene"].astype(str)
    single_gene_names = expression_data[~expression_data.Gene.str.contains("//")].reset_index(drop=True)
    multiple_gene_names = expression_data[expression_data.Gene.str.contains("//")].reset_index(drop=True)
    breaks_gene_names = pd.DataFrame(columns=["Gene", "Data"])
    
    for index, row in multiple_gene_names.iterrows():
        for genename in row["Gene"].split("///"):
            breaks_gene_names = pd.concat(
                [
                    breaks_gene_names,
                    pd.DataFrame([{"Gene": genename, "Data": row["Data"]}]),
                ],
                axis=0,
                ignore_index=True,
            )
    
    gene_expressions = pd.concat(
        [single_gene_names, breaks_gene_names], axis=0, ignore_index=True
    )
    gene_expressions.set_index("Gene", inplace=True)
    
    return gene_expressions


@contextlib.contextmanager
def suppress_stdout() -> Iterator[None]:
    with io.StringIO() as buffer:
        try:
            sys.stdout = buffer
            yield
        finally:
            sys.stdout = sys.__stdout__
