import cellxgene_census

census = cellxgene_census.open_soma(census_version="2025-11-08")

obs = cellxgene_census.get_obs(
    census,
    organism="homo_sapiens",
    value_filter='tissue_general == "brain"',
    column_names=["cell_type","tissue","tissue_general","disease"]
)

print(obs["cell_type"].value_counts().head(50))
print(obs["disease"].value_counts().head(50))

census.close()
