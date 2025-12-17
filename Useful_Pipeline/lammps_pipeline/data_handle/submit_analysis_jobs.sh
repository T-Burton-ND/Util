#!/bin/bash
set -e

echo "Submitting autocorrelations..."
qsub ./data_handle/job_autocorrelation.sh

echo "Submitting shell finder..."
qsub ./data_handle/job_find_shells.sh

echo "Submitting half-life analysis..."
qsub ./data_handle/job_half_life.sh

echo "Submitting MSD generator..."
qsub ./data_handle/job_msd_gen.sh

echo "All jobs submitted."
