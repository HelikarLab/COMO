import asyncio
from bioservices import BioDBNet
import pandas as pd
from .input_database import InputDatabase
from .output_database import OutputDatabase
from .taxon_ids import TaxonIDs


async def _async_fetch_info(
        biodbnet: BioDBNet,
        event_loop: asyncio.AbstractEventLoop,
        input_values: list[str],
        input_db: str,
        output_db: list[str],
        taxon_id: int,
        delay: int = 5,
):
    print(f"Input ")
    conversion = await event_loop.run_in_executor(
        None,  # Defaults to ThreadPoolExecutor, uses threads instead of processes. No need to modify
        biodbnet.db2db,  # The function to call
        input_db,  # The following are arguments passed to the function
        output_db,
        input_values,
        taxon_id
    )

    # If the above db2db conversion didn't work, try again until it does
    while not isinstance(conversion, pd.DataFrame):
        print(f"\nToo many requests to BioDBNet, waiting {delay} seconds and trying again.")
        await asyncio.sleep(delay)
        database_convert = await event_loop.run_in_executor(
            None,  # Defaults to ThreadPoolExecutor, uses threads instead of processes. No need to modify
            _async_fetch_info,  # The function to call
            biodbnet,
            event_loop,
            input_values,
            input_db,
            output_db,
            taxon_id,
            delay
        )

    return conversion


async def _fetch_gene_info_manager(tasks: list[asyncio.Task], batch_length: int):
    results: list[int] = []

    index: int = 0
    for task in asyncio.as_completed(tasks):
        result = await task

        # Multiply by batch length to get the approximate number of true genes being converted
        print(f"\rCollecting genes... {(index + 1) * batch_length} of {len(tasks) * batch_length} finished", end="")
        results.append(result)
        index += 1
    print()
    return results


def fetch_gene_info(
        input_values: list[str],
        input_db: InputDatabase,
        output_db: list[OutputDatabase] = None,
        taxon_id: TaxonIDs | int = TaxonIDs.HOMO_SAPIENS,
        delay: int = 5
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
    
    input_db_value = input_db.value
    if output_db is None:
        
        output_db: list = [
            OutputDatabase.GENE_SYMBOL.value,
            OutputDatabase.GENE_ID.value,
            OutputDatabase.CHROMOSOMAL_LOCATION.value
        ]
    else:
        output_db: list = [i.value for i in output_db]

    biodbnet = BioDBNet()
    dataframe_maps: pd.DataFrame = pd.DataFrame([], columns=output_db)
    dataframe_maps.index.name = input_db.value

    if type(taxon_id) == TaxonIDs:
        taxon_id_value = taxon_id.value
        if taxon_id == TaxonIDs.HOMO_SAPIENS:
            batch_len: int = 500
        else:
            batch_len: int = 300
    else:
        taxon_id_value = taxon_id
        batch_len: int = 300

    # Create a list of tasks to be awaited
    event_loop = asyncio.new_event_loop()
    asyncio.set_event_loop(event_loop)
    async_tasks = []
    for i in range(0, len(input_values), batch_len):
        # Define an upper range of values to take from input_values
        upper_range = min(i + batch_len, len(input_values))
        task = event_loop.create_task(
            _async_fetch_info(
                biodbnet=biodbnet,
                input_values=input_values[i:upper_range],
                input_db=input_db_value,
                output_db=output_db,
                taxon_id=taxon_id_value,
                delay=delay,
                event_loop=event_loop
            )
        )
    
        async_tasks.append(task)

    # database_convert = event_loop.run_until_complete(asyncio.gather(*async_tasks))
    database_convert = event_loop.run_until_complete(_fetch_gene_info_manager(tasks=async_tasks, batch_length=batch_len))
    event_loop.close()  # Close the event loop to free resources

    # Loop over database_convert to concat them into dataframe_maps
    for i, df in enumerate(database_convert):
        print(f"\rConcatenating dataframes... {i + 1} of {len(database_convert)}", end="")
        dataframe_maps = pd.concat([dataframe_maps, df], sort=False)
    print("")
    return dataframe_maps
