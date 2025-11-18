The database file would normally have a folder called 'butppg'; however, I cannot upload it. 

The .mat files are tables:

AllData has all 3888 signals.
T48 has the first 48 signals (original dataset)

new_T has 51 files used in stats.m
F1 is the result of the filter comparison made in in stats.m

no_filter is 3888 signals with no filter, but still has pre and post processing
butter_filter is 3888 signals with Butterworth filter
IIR_filter is 3888 signals filtered with IIR filtering
cheby_filter is 3888 signals filtered with Chebyshev filter
