#!/usr/bin/env python3

import sys
import os
import pymzml

from smiter.synthetic_mzml import write_scans


def main(mzml_path, file_name, start_scan, stop_scan):
    '''
    Write filtered mzML from file in the defined spectrum ID range.

    Please define full paths to input and output mzml file. Start and stop scan 
    have to be integers. If the start_scan is no MS1, the next MS1 scan will
    be used as first scan.

    usage:
        ./write_filtered_mzml.py mzml_file file_name_to_write start_scan stop_scan

    '''
    mzml_basename = mzml_path
    dirname = os.path.dirname(mzml_path)

    # [(MS1,[MS2,MS2]),...]
    scans = []

    run = pymzml.run.Reader(mzml_path)

    id_range = [n for n in range(int(start_scan),int(stop_scan),1)]
    num_scans = len(id_range)
    scan_set = None
    total_scans = 0
    for spec_pos, spec_id in enumerate(id_range):
        spectrum = run[spec_id]
        ms_level = spectrum.ms_level
        # update scan object
        spectrum.id = spectrum.ID
        spectrum.retention_time = spectrum.scan_time[0]
        print(spec_id)
        if ms_level == 1:
            if scan_set is not None:
                scans.append(scan_set)
            scan_set = [ spectrum, [] ]
            total_scans+=1
        else:
            spectrum.precursor_mz = spectrum.selected_precursors[0]['mz']
            
            precursor_i = spectrum.selected_precursors[0].get('i',0)
            spectrum.precursor_i = precursor_i
            
            charge = spectrum.selected_precursors[0].get('charge',1)
            spectrum.precursor_charge = charge
            total_scans += 1
            if scan_set is not None:
                scan_set[1].append(spectrum)
            # write rest of MS2 if no MS1 scans follows anymore
            if spec_pos + 1 == num_scans:
                scans.append(scan_set)

    write_scans(
        file_name,
        scans
    )


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(main.__doc__)
        exit()
    else:
        main(*sys.argv[1:])
