import pandas as pd
import seaborn as sns
sns.set()

# import csv
raw_df = pd.read_csv("owid-covid-data.csv", delimiter=",")

# transform date format
raw_df['date'] = pd.to_datetime(raw_df['date'])

# select country
df = raw_df.loc[(raw_df["location"] == "Switzerland")].copy()

# select time range
start_date = '2020-10-01'
end_date = '2021-02-28'
mask = (df['date'] > start_date) & (df['date'] <= end_date)
df = df.loc[mask]
