#!/bin/sh
# get data dictionary for the olink cohort which consists of individuals with olink data only 

dx extract_dataset \
        record-Gk8xFpjJgZk0BKJp0P9ZyqG0 \
        -ddd \
        --delimiter ,
