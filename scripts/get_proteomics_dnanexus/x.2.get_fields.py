# get field names for protein expression data. this is to get the whole dataset but it requires too much memory and need to run using spark cluster
import pandas as pd

data_dict_df = pd.read_csv('olink_samples.data_dictionary.csv')
# targets = pd.read_csv('ukbb_available_targets.csv') use these targets instead since it is too big otherwise

field_names = list(data_dict_df.loc[data_dict_df['entity'] == 'olink_instance_0', "name"].values)
field_names_str = [f"olink_instance_0.{f}" for f in field_names]
field_names_query = ",".join(field_names_str)

with open('samples_fields_list.txt', 'a') as the_file:
    the_file.write(field_names_query)
