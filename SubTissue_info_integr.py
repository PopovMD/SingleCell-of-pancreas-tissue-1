import pandas as pd

metadata = pd.read_csv('/Users/svetlana/Desktop/SingleCell_test-task_for_lab/metadata_FACS.csv')

pancrMD = metadata[metadata['tissue'] == 'Pancreas']

subtissues = dict(zip(pancrMD['plate.barcode'], pancrMD['subtissue']))



pancr_counts = pd.read_csv('/Users/svetlana/Desktop/SingleCell_test-task_for_lab/pancreas.raw.counts/Pancreas-counts.csv')

names = pancr_counts.columns.to_list()

for key, value in subtissues.items():
  for i in range(len(names)):
    if key in names[i]:
      names[i] = names[i] + '_' + value

pancr_counts.columns = names

pancr_counts.to_csv('pancreas_renamed_subtissue.csv', index = False)

    


