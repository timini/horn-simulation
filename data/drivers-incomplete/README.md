# Compression drivers with incomplete public data

`compression-drivers.json` lists compression drivers whose manufacturer pages or datasheets give the exit diameter, frequency range, power rating and sensitivity, but not the full Thiele-Small set (Fs, Re, Bl, Sd, Mms, Le and the Qs) that the driver model needs. Each record lists the values that were found, with the source URL and retrieval date, and names the missing fields in `missing_ts_fields`. No value was estimated.

These records are not in a driver database directory and the pipeline cannot load them. They record what is still needed: measure the missing parameters (an impedance sweep with an added-mass or closed-box step gives Fs, Re, Le, Qms, Qes, Bl and Mms; Sd needs the diaphragm dimensions) or obtain them from the manufacturer, then enter the driver with `horn-add-driver`.

Four compression drivers with complete published parameters are in `data/drivers-curated/Peerless`.
