import pandas as pd
from scipy.optimize import curve_fit

def func_4pl(x, a, b, c, d):
    return d + ((a - d) / (1 + (pow(x / c, b))))

def inverse_4pl(y, a, b, c, d):
    res = c * ((a - y) / (y - d)) ** (1 / b)
    if y <= a:
        return 0
    return res

df = pd.read_csv(input())
df2 = pd.read_csv(input())

df['AverageOD'] = (df['OD1'] + df['OD2']) / 2
average_blank = df.iloc[-1]['AverageOD']
df['CorrectedOD'] = df['AverageOD'] - average_blank


df2['AverageOD'] = (df2['OD1'] + df2['OD2']) / 2
x = df2['Concentration (ng/ml)']
y = df2['AverageOD']
params, _ = curve_fit(func_4pl, x, y)
A,B,C,D = params[0], params[1], params[2], params[3]

concentration = df['CorrectedOD'].apply(lambda y: inverse_4pl(y, A, B, C, D))
df.insert(len(df.columns), 'Concentration (ng/ml)', concentration)
detection = df['Concentration (ng/ml)'].apply(lambda x: True if x > 1 else False)
df.insert(len(df.columns), 'Protein detectable', detection)
df.to_csv('elisa_result.csv')