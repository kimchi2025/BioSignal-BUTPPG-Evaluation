# Done in python because it was easier
import pandas as pd

# loading dataframes
hr = pd.read_csv('quality-hr-ann.csv')
subj = pd.read_csv('subject-info.csv')

# extracting ID, ear/finger, and motion columns from subject dataframe then merging with hr dataframe
merged = subj[['ID', 'Ear/finger', 'Motion']].merge(hr, on='ID')

# saving merged dataframe as a new csv file to use in MATLAB
merged.to_csv('subj_hr_info.csv', index=False)
