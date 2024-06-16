#!/bin/bash


FIELDS=$(cat samples_fields_list.txt)
# FIELDS=$(cat ukbb_available_targets_query.txt) #this is to only get genes of interest

dx extract_dataset \
        record-Gk8xFpjJgZk0BKJp0P9ZyqG0 \
        --fields $FIELDS \
        --delimiter , \
        --output get_expressions_query.txt \
	--sql


