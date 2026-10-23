# Reference-genome validation

`expected.tsv` lists public reference genomes (NCBI accessions) and their
expected ST and tox status. The *Validate* workflow downloads them from NCBI,
types them with the current database, and fails if any result differs. It
runs automatically when the database or typing script changes, so every
database update is checked against known genomes.

To add a genome, append a row: `sample  accession  expected_ST  positive|negative`.
Only add publicly available genomes - never patient isolates.
