import smiter
import matplotlib.pyplot as plt


def main():
    peak_props = {
        "inosine": {
            "trivial_name": "inosine",
            "chemical_formula": "+C(10)H(12)N(4)O(5)",
            "charge": 2,
            "scan_start_time": .0,
            "peak_width": 2,  # minutes
            "peak_function": "gauss_tail",
            "peak_scaling_factor": 1e3,
        }
    }
    mzml_params = {
        'gradient_length': 3,
        'min_intensity': 50,
        'isolation_window_width': 0.7,
        "ms_rt_diff": 0.001,  # in minutes
    }

    noise_injector = smiter.noise_functions.JamssNoiseInjector()
    fragmentor = smiter.fragmentation_functions.NucleosideFragmentor()

    interval_tree = smiter.synthetic_mzml.genereate_interval_tree(peak_props)
    trivial_names = {
        val["chemical_formula"]: key for key, val in peak_props.items()
    }
    isotopologue_lib = smiter.synthetic_mzml.generate_molecule_isotopologue_lib(
        peak_props, trivial_names=trivial_names
    )

    fig, (ax1, ax2) = plt.subplots(1, 2)
    scans_noise, scan_dict = smiter.synthetic_mzml.generate_scans(
        isotopologue_lib,
        peak_props,
        interval_tree,
        fragmentor,
        noise_injector,
        mzml_params,
    )
    t = [scan['rt'] for scan, _ in scans_noise]
    i = [sum(scan['i']) for scan, _ in scans_noise]
    ax1.plot(t, i)
    ax1.set_title("With noise")
    # plt.show()

    noise_injector = smiter.noise_functions.UniformNoiseInjector(
        dropout=0,
        ppm_noise=0,
        intensity_noise=0,
    )
    scans, scan_dict = smiter.synthetic_mzml.generate_scans(
        isotopologue_lib,
        peak_props,
        interval_tree,
        fragmentor,
        noise_injector,
        mzml_params,
    )
    t = [scan['rt'] for scan, _ in scans]
    i = [sum(scan['i']) for scan, _ in scans]
    ax2.plot(t, i)
    ax2.set_title("Without noise")
    plt.show()

    for s, _ in scans_noise:
        plt.title(f'Spectrum {s.id}')
        plt.bar(s.mz, s.i)
        plt.show()




if __name__ == '__main__':
    main()
