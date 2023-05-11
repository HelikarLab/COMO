import asyncio
import pandas as pd
from bioservices import BioDBNet
import bioservices

try:
    from input_database import InputDatabase
    from output_database import OutputDatabase
    from taxon_ids import TaxonIDs
    from utilities import suppress_stdout
except ImportError:
    from .input_database import InputDatabase
    from .output_database import OutputDatabase
    from .taxon_ids import TaxonIDs
    
    import sys
    sys.path.append("..")
    from utilities import suppress_stdout

async def _async_fetch_info(
        biodbnet: BioDBNet,
        event_loop: asyncio.AbstractEventLoop,
        semaphore: asyncio.Semaphore,
        input_values: list[str],
        input_db: str,
        output_db: list[str],
        taxon_id: int,
        delay: int = 10
) -> pd.DataFrame:
    await semaphore.acquire()
    conversion = await asyncio.to_thread(
        biodbnet.db2db,
        input_db,
        output_db,
        input_values,
        taxon_id
    )
    semaphore.release()
    
    # If the above db2db conversion didn't work, try again until it does
    if not isinstance(conversion, pd.DataFrame):
        # Errors will occur on a timeouut. If this happens, split our working dataset in two and try again
        first_set: list[str] = input_values[:len(input_values) // 2]
        second_set: list[str] = input_values[len(input_values) // 2:]
        
        await asyncio.sleep(delay)
        first_conversion: pd.DataFrame = await _async_fetch_info(
            biodbnet, event_loop, semaphore,
            first_set, input_db, output_db,
            taxon_id, delay
        )
        second_conversion: pd.DataFrame = await _async_fetch_info(
            biodbnet, event_loop, semaphore,
            second_set, input_db, output_db,
            taxon_id, delay,
        )
        return pd.concat([first_conversion, second_conversion])
        
    return conversion


async def _fetch_gene_info_manager(tasks: list[asyncio.Task[pd.DataFrame]], batch_length: int) -> list[pd.DataFrame]:
    results: list[pd.DataFrame] = []

    print("Collecting genes... ", end="")
    
    task: asyncio.Future[pd.DataFrame]
    for i, task in enumerate(asyncio.as_completed(tasks)):
        results.append(await task)
        print(f"\rCollecting genes... {(i + 1) * batch_length} of ~{len(tasks) * batch_length} finished", end="")

    return results

def fetch_gene_info(
        input_values: list[str],
        input_db: InputDatabase,
        output_db: OutputDatabase | list[OutputDatabase] = None,
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
    batch_length: int = 100
    output_db_values: list[str]
    if output_db is None:
        output_db_values = [
            OutputDatabase.GENE_SYMBOL.value,
            OutputDatabase.GENE_ID.value,
            OutputDatabase.CHROMOSOMAL_LOCATION.value
        ]
    elif isinstance(output_db, OutputDatabase):
        output_db_values = [output_db.value]
    else:
        output_db_values = [str(i.value) for i in output_db]

    if type(taxon_id) == TaxonIDs:
        taxon_id_value = taxon_id.value
    else:
        taxon_id_value = taxon_id

    using_cache: bool = True
    biodbnet: BioDBNet
    try:
        biodbnet = BioDBNet(cache=True, verbose=False)
        _ = biodbnet.services.session.settings.cache_control
    except ImportError:  # Error in case sqlite is not found for some reason
        biodbnet = BioDBNet(cache=False, verbose=False)
        using_cache = False
    biodbnet.services.TIMEOUT = 60
    
    if using_cache:
        print(f"Using cache for BioDBNet")
    else:
        print(f"Unable to set cache for BioDBNet")

    dataframe_maps: pd.DataFrame = pd.DataFrame([], columns=output_db_values)
    dataframe_maps.index.name = input_db.value
    
    # Create a list of tasks to be awaited
    event_loop = asyncio.new_event_loop()
    asyncio.set_event_loop(event_loop)
    async_tasks = []
    semaphore: asyncio.Semaphore = asyncio.Semaphore(15)
    for i in range(0, len(input_values), batch_length):
        # Define an upper range of values to take from input_values
        upper_range = min(i + batch_length, len(input_values))
        task = event_loop.create_task(
            _async_fetch_info(
                biodbnet=biodbnet,
                semaphore=semaphore,
                input_values=input_values[i:upper_range],
                input_db=input_db_value,
                output_db=output_db_values,
                taxon_id=taxon_id_value,
                delay=delay,
                event_loop=event_loop
            )
        )
    
        async_tasks.append(task)

    database_convert = event_loop.run_until_complete(_fetch_gene_info_manager(tasks=async_tasks, batch_length=batch_length))
    event_loop.close()  # Close the event loop to free resources

    # Loop over database_convert to concat them into dataframe_maps
    print("")
    for i, df in enumerate(database_convert):
        print(f"Concatenating dataframes... {i + 1} of {len(database_convert)}" + " " * 50, end="\r")
        dataframe_maps = pd.concat([dataframe_maps, df], sort=False)
    print("")
    return dataframe_maps
