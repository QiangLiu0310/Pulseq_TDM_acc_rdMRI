
1. Edit `set_experimental_params.m` 
   1. pf_ky, readout .mod file, file names, etc

2. Write raw data to .h5 files
   ```
   datfiles_to_h5;
   ```

3. Get ghost correction parameters `a`
   ```
   get_ghost_calibration_data;
   ```
   Change `del` in `set_experimental_params.m` as needed to avoid phase wraps in ghost correction step

4. Get 2d reference data
   ```
   get_acs_data;
   ```
4. slice GRAPPA recon
   1. Set desired frame(s) to recon in `recon_timeseries.m`
   Then:
   ```
   recon_timeseries;
   ```





