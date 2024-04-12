import pandas as pd
import sys

adult_list = [
    'early adult stage',
    '29-year-old human stage',
    '42-year-old human stage',
    '50-year-old human stage',
    '60-year-old human stage',
    '38-year-old human stage',
    '46-year-old human stage',
    '54-year-old human stage',
    '56-year-old human stage',
    '61-year-old human stage',
    '62-year-old human stage',
    '63-year-old human stage',
    '67-year-old human stage',
    '68-year-old human stage',
    '70-year-old human stage',
    '73-year-old human stage',
    '75-year-old human stage',
    '78-year-old human stage',
    '82-year-old human stage',
    '24-year-old human stage',
    '35-year-old human stage',
    '36-year-old human stage',
    '41-year-old human stage',
    '49-year-old human stage',
    '51-year-old human stage',
    '53-year-old human stage',
    '57-year-old human stage',
    '58-year-old human stage',
    '59-year-old human stage',
    '71-year-old human stage',
    '45-year-old human stage',
    '21-year-old human stage',
    '22-year-old human stage',
    '25-44 year-old human stage',
    '26-year-old human stage',
    '28-year-old human stage',
    '30-year-old human stage',
    '31-year-old human stage',
    '32-year-old human stage',
    '34-year-old human stage',
    '37-year-old human stage',
    '39-year-old human stage',
    '44-year-old human stage',
    'human adult stage',
    '52-year-old human stage',
    'human late adulthood stage',
    '64-year-old human stage',
    '65-year-old human stage',
    '72-year-old human stage',
    '55-year-old human stage',
    '66-year-old human stage',
    '74-year-old human stage',
    '40-year-old human stage',
    '69-year-old human stage',
    '76-year-old human stage',
    '79-year-old human stage',
    '80 year-old and over human stage',
    '81-year-old human stage',
    '87-year-old human stage',
    '89-year-old human stage',
    'prime adult stage',
    'late adult stage',
    '48-year-old human stage',
    'fifth decade human stage',
    'fourth decade human stage',
    'seventh decade human stage',
    'sixth decade human stage',
    'third decade human stage',
    '47-year-old human stage',
    '77-year-old human stage',
    'eighth decade human stage',
    '43-year-old human stage',
    '19-year-old human stage',
    '20-year-old human stage',
    '25-year-old human stage',
    '33-year-old human stage',
    '80-year-old human stage',
    '83-year-old human stage',
    '84-year-old human stage',
    '85-year-old human stage',
    '86-year-old human stage',
    '88-year-old human stage',
    '90-year-old human stage',
    '91-year-old human stage',
    '93-year-old human stage',
    '94-year-old human stage',
    '95-year-old human stage',
    '96-year-old human stage',
    '97-year-old human stage',
    '65-79 year-old human stage',
    'ninth decade human stage',
    'human early adulthood stage',
    'human middle aged stage',
    'human aged stage',
    '18-year-old human stage'
]

file_path = sys.argv[1]
df = pd.read_csv(file_path)

unique = list(df['development_stage'].unique())
d1 = list(set(adult_list) - set(unique))
d2= list(set(unique) - set(adult_list))
print('items in adult_list not in unique:', d1)
print('items in unique not in adult_list:', d2)

print(d2)
fil = pd.DataFrame()
fil['development_stage'] = d2

young = df.merge(fil, on = 'development_stage')['dataset_id'].unique()
print(young)
print(len(young))

adults = df[~df['dataset_id'].isin(young)]

print(adults)
adults_uni = adults.drop_duplicates(subset='dataset_id',keep='first')
adults_uni[['dataset_id']].to_csv('adult_dataset_ids.csv', index=False)

print('reading in giant dataframe')
big_df = pd.read_csv('../../expression-summary-full-03-11-24.csv')

print('read')
adults_unid = adults_uni[['dataset_id']]

print('removing young IDs')
big_df = bif_df[~big_df['dataset_id'].isin(young)]

print('writing file out')
big_df.to_csv('../../expression-summary-adults-04-11-24.csv', index=False)
