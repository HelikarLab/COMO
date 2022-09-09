import asyncio
from bioservices import BioDBNet
import pandas as pd
from .input_database import InputDatabase
from .output_database import OutputDatabase
from.taxon_ids import TaxonIDs


async def _async_fetch_info(
        biodbnet: BioDBNet,
        input_values: list[str],
        input_db: str,
        output_db: list[str],
        taxon_id: TaxonIDs | int,
        delay: int
):

    if isinstance(taxon_id, TaxonIDs):
        taxon_id_value = taxon_id.value
    else:
        taxon_id_value = taxon_id

    database_convert = await asyncio.to_thread(
        biodbnet.db2db,  # The function to call
        input_db,       # The following are arguments passed to the function
        output_db,
        input_values,
        taxon_id_value
    )

    # If the above db2db conversion didn't work, try again until it does
    while not isinstance(database_convert, pd.DataFrame):
        print(f"\nToo many requests to BioDBNet, waiting {delay} seconds and trying again.")
        await asyncio.sleep(delay)
        database_convert = await asyncio.to_thread(
            _async_fetch_info,  # The function to call
            biodbnet,
            input_values,
            input_db,
            output_db,
            taxon_id,
            delay
        )

    return database_convert


async def _async_fetch_info_wrapper(coroutines: list):
    return await asyncio.gather(*coroutines)


def fetch_gene_info(
        input_values: list[str],
        input_db: InputDatabase,
        output_db: list[OutputDatabase] | list[str],
        taxon_id: TaxonIDs | int = TaxonIDs.HUMAN,
        delay: int = 10
) -> pd.DataFrame:
    """
    This function returns a dataframe with important gene information for future operations in MADRID.

    Fetch gene information from BioDBNet
    :param input_values: A list of genes in "input_db" format
    :param input_db: The input database to use (default: "Ensembl Gene ID")
    :param output_db: The output format to use (default: ["Gene Symbol", "Gene ID", "Chromosomal Location"])
    :param delay: The delay in seconds to wait before trying again if bioDBnet is busy (default: 15)
    :param taxon_id: The taxon ID to use (default: 9606)

    :return: A dataframe with specified columns as "output_db" (Default is HUGO symbol, Entrez ID, and chromosome start and end positions)
    """
    input_db = input_db.value
    if output_db is None:
        output_db: list[str] = [
            OutputDatabase.GENE_SYMBOL.value,
            OutputDatabase.GENE_ID.value,
            OutputDatabase.CHROMOSOMAL_LOCATION.value
        ]
    else:
        output_db: list[str] = [i.value for i in output_db]

    biodbnet = BioDBNet()
    dataframe_maps: pd.DataFrame = pd.DataFrame([], columns=output_db)
    dataframe_maps.index.name = input_db

    if taxon_id == TaxonIDs.HUMAN:
        batch_len = 500
    else:
        batch_len = 300

    # Create a list of tasks to be awaited
    async_coroutines = []
    for i in range(0, len(input_values), batch_len):
        # Define an upper range of values to take from input_values
        upper_range = min(i + batch_len, len(input_values))

        async_coroutines.append(
            _async_fetch_info(
                biodbnet=biodbnet,
                input_values=input_values[i:upper_range],
                input_db=input_db,
                output_db=output_db,
                taxon_id=taxon_id,
                delay=delay
            )
        )

    # Await all the tasks, must use [0] to get completed results. [1] is a set of pending tasks (i.e., empty).
    database_convert = asyncio.run(_async_fetch_info_wrapper(async_coroutines))

    # Loop over database_convert to concat them into dataframe_maps
    for df in database_convert:
        dataframe_maps = pd.concat([dataframe_maps, df], sort=False)

    return dataframe_maps
