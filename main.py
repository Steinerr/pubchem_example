import math
import pubchempy as pcp
import matplotlib.pyplot as plt
import pandas as pd

cids_95 = pcp.get_cids(
    'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
    'smiles',
    searchtype='similarity',
    Threshold=95
)
print(f'len cids for 95: {len(cids_95)}')

cids_80 = pcp.get_cids(
    'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
    'smiles',
    searchtype='similarity',
    Threshold=80
)
print(f'len cids for 80: {len(cids_80)}')

cids_70 = pcp.get_cids(
    'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C',
    'smiles',
    searchtype='similarity',
    Threshold=70
)
print(f'len cids for 70: {len(cids_70)}')

# график зависимости количества результатов поиска от порога схожести
threshholds = [95, 80, 70]
threshholds_lengths = [len(cids_95), len(cids_80), len(cids_70)]

plt.plot(threshholds, threshholds_lengths)
plt.xlabel('Threshhold similarity')
plt.ylabel('Count of results')
plt.savefig('threshholds_length.png')

# pandas dataframe c четырьмя колонками
cids_95_df = pcp.get_properties(
  ('XlogP','TPSA','MolecularWeight','Volume3D'), 
  cids_95, 
  namespace='cid', 
  as_dataframe=True
)

cids_80_df = pcp.get_properties(
  ('XlogP','TPSA','MolecularWeight','Volume3D'), 
  cids_80, 
  namespace='cid', 
  as_dataframe=True
)

listkey_count = 5000
pages = math.ceil(len(cids_70) / listkey_count)
print(f'Found {pages} pages')
cids_70_dfs = []
for page_number in range(1, pages + 1):
  print(f'Downloading page #{page_number} from {pages}')
  df = pcp.get_properties(
    ('XlogP','TPSA','MolecularWeight','Volume3D'), 
    cids_70[(page_number - 1) * listkey_count:page_number * listkey_count], 
    namespace='cid', 
    as_dataframe=True
  )
  cids_70_dfs.append(df)

cids_70_df = pd.concat(cids_70_dfs)

# pandas dataframe без строк имеющих значения None
print('Dropping NaNs')
cids_95_df = cids_95_df.dropna(axis=0, how='any')
cids_80_df = cids_80_df.dropna(axis=0, how='any')
cids_70_df = cids_70_df.dropna(axis=0, how='any')

# гистограмма
ax = cids_95_df.plot.hist(bins=100, alpha=0.5)
fig = ax.get_figure()
fig.savefig('hist_95_full.png')

ax = cids_80_df.plot.hist(bins=100, alpha=0.5)
fig = ax.get_figure()
fig.savefig('hist_80_full.png')

ax = cids_70_df.plot.hist(bins=100, alpha=0.5)
fig = ax.get_figure()
fig.savefig('hist_70_full.png')

# Средние значения по колонкам
print('Averages:')
print(f'Threshhold 95: {cids_95_df.mean(axis=0)}')
print(f'Threshhold 80: {cids_80_df.mean(axis=0)}')
print(f'Threshhold 70: {cids_70_df.mean(axis=0)}')

# Медиана
print('Medians:')
print(f'Threshhold 95: {cids_95_df.median(axis=0)}')
print(f'Threshhold 80: {cids_80_df.median(axis=0)}')
print(f'Threshhold 70: {cids_70_df.median(axis=0)}')

# Матрицы рассеяния
c = pd.plotting.scatter_matrix(cids_95_df, alpha=0.2)
plt.savefig('matrix_95.png')

c = pd.plotting.scatter_matrix(cids_80_df, alpha=0.2)
plt.savefig('matrix_80.png')

c = pd.plotting.scatter_matrix(cids_70_df, alpha=0.2)
plt.savefig('matrix_70.png')

# Корреляция по Пирсону
print(f'Corr Threshold 95: {cids_95_df.corr(method = "pearson")}')
print(f'Corr Threshold 80: {cids_80_df.corr(method = "pearson")}')
print(f'Corr Threshold 70: {cids_70_df.corr(method = "pearson")}')
