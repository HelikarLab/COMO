from bioservices import BioDBNet

uniprot: list[str] = (
    open(
        "/Users/joshl/Library/Application Support/JetBrains/PyCharm2022.1/scratches/uniprot_ids.txt"
    )
    .read()
    .splitlines()
)

biodbnet = BioDBNet()

entrez_id = biodbnet.db2db(
    input_db="UniProt Accession",
    output_db="Gene ID",
    input_values=uniprot,
)

for i in entrez_id["Gene ID"]:
    print(i)
