# Subject 1 scan, UM UHP Jan 24, 2024


Other recent changes:
  * getsmspulse.m: put first slice at -mb/2 (remove shift of +sliceSep/2)
     * writeEPI.m: accordingly, remove shift of -sliceSep/2 (rf.freqOffset)
  * set maxView = np x etl when mb>1, etl otherwise
  * Fixed slice offset for mb=1
  * etl=72 (multiple of mb=6)
  * TR=800ms
  * Add RF spoiling
  * add fat sat as default
  * Interleaved partition ordering, with last two even shots swapped
  * Remove arg.spoilersOn option (always on)

