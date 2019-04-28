import pandas as pd


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


inp = [{'start': 60, 'end': 70, 'x': 1},
       {'start': 90, 'end': 210, 'x': 2},
       {'start': 80, 'end': 110, 'x': 2},
       {'start': 110, 'end': 190, 'x': 2},
       {'start': 190, 'end': 210, 'x': 2},
       {'start': 220, 'end': 250, 'x': 2}]
df = pd.DataFrame(inp)
print(df.head())

start_idx = 100
end_idx = 200

df_overlap = df[((df['end'] >= start_idx) & (df['end'] <= end_idx))
                | ((df['start'] >= start_idx) & (df['start'] <= end_idx))
                | ((df['start'] <= start_idx) & (df['end'] >= end_idx))]

print(df_overlap)


x_average = 0
total_covered = 0
for index, row in df_overlap.iterrows():
    overlap_window = getOverlap([start_idx, end_idx], [row['start'], row['end']])
    x_average += row['x']*overlap_window/(row['end'] - row['start'])
    total_covered += overlap_window

print(x_average/total_covered)


check = ( ((2/120)*100) + ((2/30)*10) + ((2/80)*80)+ ((2/20)*10) ) /(100 + 10 + 80 + 10)
print(check)




