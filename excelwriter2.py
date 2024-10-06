import pandas as pd
from pandas import ExcelWriter

writer = ExcelWriter('pandas.xlsx', mode='w', engine='xlsxwriter')
worksheet = writer.sheets['Sheet1']
workbook = writer.book
header_format = workbook.add_format({
  'bold': True,
  'text_wrap': True,
  'valign': 'top',
  'fg_color': '#D7E4BB',
  'border': 1
})
# write the header separately so you can format it
for col_num, value in enumerate(['Manager', 'Policy Name', 'Numbers', 'Source IP', 'Destination IP']):
  worksheet.write(0, col_num + 1, value, header_format)

df = pd.read_excel(file)
# if you only want specific columns
df = df[['Manager', 'Policy Name', 'Numbers', 'Source', 'Destination']]
df.to_excel(writer, sheet_name="Holi", index=False, header=False, startrow=1)
writer.save()