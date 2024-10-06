import pandas as pd

df1 = pd.DataFrame({'Data': ['a', 'b', 'c', 'd']})

df2 = pd.DataFrame({'Data': [1, 2, 3, 4]})

df3 = pd.DataFrame({'Data': [1.1, 1.2, 1.3, 1.4]})

writer = pd.ExcelWriter('multiple.xlsx', engine='xlsxwriter')

df1.to_excel(writer, sheet_name='Sheeta')
df2.to_excel(writer, sheet_name='Sheetb')
df3.to_excel(writer, sheet_name='Sheetc')

writer.save()
