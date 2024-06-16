#!/bin/bash
# this scripts extract data fields that are not present in our current basket.

dx extract_dataset \
	record-Gk420f8J24QGxpj5yqYZgjkG \
	--fields participant.eid,participant.p23105_i0,participant.p399_i0_a1,participant.p399_i0_a2,participant.p399_i0_a3,participant.p10137_i0_a1,participant.p10137_i0_a2,participant.p10137_i0_a3,participant.p4282_i0,participant.p20016_i0 \
	--delimiter , \
	--output additional_fields.csv
	
