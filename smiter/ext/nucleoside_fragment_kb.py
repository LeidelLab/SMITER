"""Mapping nucleoside names to its fragments."""

KB_FRAGMENTATION_INFO: dict = {
    "uridine": {
        # U
        "comments": "precursor can not be seen in MS2 in some cases",
        "references": ["Jora et al. 2018"],
        "fragments": {
            "uridine -Ribose": {"formula": "C(4)H(4)N(2)O(2)"},  # 113.0345538436
            "uridine -Ribose -N(1)H(3)": {
                "formula": "C(4)H(1)N(1)O(2)",  # 96.0080047430
                "hcd": True,
            },
        },
        "precursors": {
            "uridine": {},
            "uridine dimer": {},
        },
    },
    "5-methyluridine": {
        # m5U
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "5-methyluridine": {},
            "5-methyluridine +N(1)H(2)": {"formula": "C(10)H(16)N(3)O(6)"},
            "5-methyluridine dimer": {},
        },
        "fragments": {
            "5-methyluridine -Ribose": {"formula": "C(5)H(6)N(2)O(2)"},
            "5-methyluridine -Ribose -N(1)H(3)": {
                "formula": "C(5)H(3)N(1)O(2)",  # 110.0236548074
                "hcd": True,
            },
            "5-methyluridine -Ribose -H(2)O(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            # '5-methyluridine -Ribose -C(1)H(1)N(1)O(1)' : {
            #     'formula': 'C(4)H(5)N(1)O(1)', #84.04
            #     'hcd':True
            # },
            # '5-methyluridine -Ribose -C(1)H(3)N(1)O(1)' : {
            #     'formula': 'C(4)H(3)N(1)O(1)', #82.028
            #     'hcd':True
            # },
            # '5-methyluridine -Ribose -C(2)H(2)N(2)O(1)' : {
            #     'formula': 'C(3)H(4)O(1)', #57.033
            #     'hcd':True
            # },
            # '5-methyluridine -Ribose -C(2)H(1)N(1)O(2)' : {
            #     'formula': 'C(3)H(5)N(1)', #56.050
            #     'hcd':True
            # },
            # '5-methyluridine -Ribose -C(2)H(3)N(1)O(2)' : {
            #     'formula': 'C(3)H(3)N(1)', #53.034
            #     'hcd':True
            # },
        },
        "exclusion_fragments": {
            "5-methyluridine -Methylated ribose": {
                "formula": "C(4)H(4)N(2)O(2)"  # 113.0345538436
            }
        },
    },
    "3-methyluridine": {
        # m3U
        "comments": "Exclusion fragments are partly m5U specific fragments",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "3-methyluridine": {},
            "3-methyluridine +N(1)H(2)": {"formula": "C(10)H(16)N(3)O(6)"},
            "3-methyluridine dimer": {},
        },
        "fragments": {
            "3-methyluridine -Ribose": {"formula": "C(5)H(6)N(2)O(2)"},
            "3-methyluridine -Ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(1)N(1)O(2)",  # 96
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "3-methyluridine -Methylated ribose": {
                "formula": "C(4)H(4)N(2)O(2)"  # 113.0345538436
            },
            "3-methyluridine -Ribose -N(1)H(3)": {
                "formula": "C(5)H(3)N(1)O(2)",  # 110.02
                "hcd": True,
            },
            "3-methyluridine -Ribose -H(2)O(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.03
                "hcd": True,
            },
        },
    },
    "2′-O-methyluridine": {
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"2′-O-methyluridine": {}},
        "fragments": {
            "2′-O-methyluridine -Methylated ribose": {
                "formula": "C(4)H(4)N(2)O(2)"  # 113.035
            },
            "2′-O-methyluridine -Methylated ribose -H(3)N(1)": {
                "formula": "C(4)H(1)N(1)O(2)",  # 96.008
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "2′-O-methyluridine -Ribose": {"formula": "C(5)H(6)N(2)O(2)"}  # 127
        },
    },
    "2′-O-methylpseudouridine": {
        "comments": "Verification required",
        "references": [],
        "precursors": {"2′-O-methylpseudouridine": {}},
        "fragments": {
            "2′-O-methylpseudouridine -H(4)O(4)": {"formula": "C(10)H(10)N(2)O(4)"}
        },
    },
    "pseudouridine": {
        "references": ["Dudley et al. 2005", "Pomerantz et al. 2005"],
        "comments": "peak 209.0556832124: ribose is doubly dehydrated; peak 191.0451185280: ribose is triply dehydrated",
        "fragments": {
            "pseudouridine -H(4)O(2)": {
                "formula": "C(9)H(8)N(2)O(4)",  # 209.0556832124
                "hcd": False,
            },
            "pseudouridine -H(6)O(3)": {
                "formula": "C(9)H(6)N(2)O(3)",  # 191.0451185280
                "hcd": False,
            },
            "pseudouridine -C(3)H(6)O(3)": {
                "formula": "C(6)H(6)N(2)O(3)",  # 155.0451185280
                "hcd": False,
            },
            "pseudouridine -C(4)H(8)O(4)": {
                "formula": "C(5)H(4)N(2)O(2)",  # 125.0345538436
                "hcd": True,
            },
            "pseudouridine -C(5)H(7)N(1)O(4)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "pseudouridine -C(5)H(9)N(1)O(5)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
            "pseudouridine -C(5)H(6)N(1)O(6)": {
                "formula": "C(4)H(5)N(1)",  # 68.0494756318
                "hcd": True,
            },
        },
        "precursors": {
            "pseudouridine K adduct": {"formula": "C(9)H(11)N(2)O(6)K(1)"},
            "pseudouridine": {},
            "pseudouridine +N(1)H(2)": {"formula": "C(9)H(14)N(3)O(6"},
            # 'pseudouridine dimer' : {},
            # 'pseudouridine K adduct' : {
            # 'formula' : 'C(9)H(11)N(2)O(6)K(1)'
            # },
            # 'pseudouridine +Ammonia' : {
            # 'formula' : 'C(9)H(15)N(3)O(6)'
            # },
        },
    },
    "1-methylpseudouridine": {
        # m1Y
        "references": ["Dudley et al. 2005"],
        "comments": "Peaks similar to pseudouridine with addtional +C(1)H(2)",
        "fragments": {
            "1-methylpseudouridine -H(4)O(2)": {
                "formula": "C(10)H(10)N(2)O(4)",  # 223.0713332768
                "hcd": False,
            },
            "1-methylpseudouridine -H(6)O(3)": {
                "formula": "C(10)H(8)N(2)O(3)",  # 205.0607685924
                "hcd": False,
            },
            "1-methylpseudouridine -C(3)H(6)O(3)": {
                "formula": "C(7)H(8)N(2)O(3)",  # 169.0607685924
                "hcd": False,
            },
            "1-methylpseudouridine -C(4)H(8)O(4)": {
                "formula": "C(6)H(6)N(2)O(2)",  # 139.0502039080
                "hcd": True,
            },
            "1-methylpseudouridine -C(5)H(9)N(1)O(5)": {
                "formula": "C(5)H(5)N(1)O(1)",  # 96.0443902518
                "hcd": True,
            },
        },
        "precursors": {
            "1-methylpseudouridine K adduct": {},
            "1-methylpseudouridine": {},
            "1-methylpseudouridine dimer": {},
        },
    },
    "3-methylpseudouridine": {
        # m3Y
        "references": ["Dudley et al. 2005"],
        "comments": "Verification required; Peaks similar to pseudouridine with addtional +C(1)H(2)",
        "fragments": {
            "3-methylpseudouridine -H(4)O(2)": {
                "formula": "C(10)H(10)N(2)O(4)",  # 223.0713332768
                "hcd": False,
            },
            "3-methylpseudouridine -H(6)O(3)": {
                "formula": "C(10)H(8)N(2)O(3)",  # 205.0607685924
                "hcd": False,
            },
            "3-methylpseudouridine -C(3)H(6)O(3)": {
                "formula": "C(7)H(8)N(2)O(3)",  # 169.0607685924
                "hcd": False,
            },
            "3-methylpseudouridine -C(4)H(8)O(4)": {
                "formula": "C(6)H(6)N(2)O(2)",  # 139.0502039080
                "hcd": True,
            },
            "3-methylpseudouridine -C(5)H(9)N(1)O(5)": {
                "formula": "C(5)H(5)N(1)O(1)",  # 96.0443902518
                "hcd": True,
            },
        },
        "precursors": {
            "3-methylpseudouridine K adduct": {},
            "3-methylpseudouridine": {},
            "3-methylpseudouridine dimer": {},
        },
    },
    "inosine": {
        # I
        "references": [],
        "comments": "Precursor ion can not be seen in MS2 in most cases",
        "fragments": {
            "inosine -Ribose": {"formula": "C(5)H(4)N(4)O(1)"}  # 137.0457872316
        },
        "precursors": {"inosine": {}, "inosine dimer": {}},
    },
    "1-methylinosine": {
        # m1I
        "references": [],
        "comments": "",
        "precursors": {"1-methylinosine": {}, "1-methylinosine dimer": {}},
        "fragments": {
            "1-methylinosine -Ribose": {
                "formula": "C(6)H(6)N(4)O(1)"  # 151.06143729597002
            }
        },
        "exclusion_fragments": {
            "1-methylinosine -Methylated ribose": {
                "formula": "C(5)H(4)N(4)O(1)"  # 137.04578723157002
            }
        },
    },
    "2′-O-methylinosine": {
        # Im
        "references": [],
        "comments": "",
        "precursors": {"2′-O-methylinosine": {}, "2′-O-methylinosine dimer": {}},
        "fragments": {
            "2′-O-methylinosine -Methylated ribose": {
                "formula": "C(5)H(4)N(4)O(1)"  # 137.04578723157002
            }
        },
        "exclusion_fragments": {
            "2′-O-methylinosine -Ribose": {
                "formula": "C(6)H(6)N(4)O(1)"  # 151.06143729597002
            }
        },
    },
    "N6-isopentenyladenosine": {
        # i6A
        "references": ["modomics.org"],
        "comments": "",
        "precursors": {"N6-isopentenyladenosine": {}},
        "fragments": {
            "N6-isopentenyladenosine -Ribose": {
                "formula": "C(10)H(13)N(5)",  # 204.1243719054
                "hcd": False,
            },
            "N6-isopentenyladenosine -Ribose -C(4)H(8)": {
                "formula": "C(6)H(5)N(5)",  # 148.0617716478
                "hcd": True,
            },
            "N6-isopentenyladenosine -Ribose -C(5)H(8)": {
                "formula": "C(5)H(5)N(5)",  # 136.0617716478
                "hcd": True,
            },
        },
    },
    "N6-formyladenosine": {
        # f6A
        "references": ["Fu et al. 2013"],
        "comments": "Verification by standard required",
        "precursors": {"N6-formyladenosine": {}},
        "fragments": {
            "N6-formyladenosine -Ribose": {
                "formula": "C(6)H(5)N(5)O(1)",  # 164.05668626777
                "hcd": False,
            },
            "N6-formyladenosine -Ribose -C(1)O(1)": {
                "formula": "C(5)H(5)N(5)",  # 136.0617716478
                "hcd": True,
            },
        },
    },
    "5-methoxyuridine": {
        # mo5U
        "comment": "Can be a co-eluting fragment (in-source?) of cm5U but with different fragments",
        "references": ["Nelson and McCloskey 1994"],
        "fragments": {
            "5-methoxyuridine -Ribose": {
                "formula": "C(5)H(6)N(2)O(3)",  # 143.0451185280
                "hcd": False,
            },
            # "5-methoxyuridine -Ribose -N(1)H(3)": {
            #     "formula": "C(5)H(3)N(1)O(3)",  # 126.0185694274
            #     "hcd": True,
            # },
        },
        "precursors": {
            "5-methoxyuridine": {},
            "5-methoxyuridine dimer": {},
            "5-methoxyuridine K adduct": {},
        },
    },
    "N6-threonylcarbamoyladenosine": {
        # t6A
        "comments": "120, 136, 162 and 281 are major peaks",
        "references": ["Nagao et al. 2017"],
        "fragments": {
            # 'N6-threonylcarbamoyladenosine -Ribose' : {
            #     'formula' : 'C(10)H(12)N(6)O(4)' # 281.0992793572
            # },
            "N6-threonylcarbamoyladenosine -Ribose -C(4)H(9)N(1)O(3)": {
                "formula": "C(6)H(3)N(5)O(1)"  # 162.0410362034
            },
            "N6-threonylcarbamoyladenosine -Ribose -C(5)H(7)N(1)O(4)": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            },
            "N6-threonylcarbamoyladenosine -Ribose -C(6)H(3)N(5)O(1)": {
                "formula": "C(4)H(9)N(1)O(3)",  # 120.0655196206
                "hcd": False,
            },
        },
        "precursors": {
            "N6-threonylcarbamoyladenosine": {
                "formula": "C(15)H(20)N(6)O(8)"  # 413.1415380948
            }
        },
    },
    "cyclic N6-threonylcarbamoyladenosine": {
        "references": ["Myauchi et al. 2013"],
        "comments": "",
        "precursors": {"cyclic N6-threonylcarbamoyladenosine": {}},
        "fragments": {
            "cyclic N6-threonylcarbamoyladenosine -Ribose": {
                "formula": "C(10)H(10)N(6)O(3)"  # 263.0887146728
            },
            "cyclic N6-threonylcarbamoyladenosine -Ribose -C(2)H(4)O(1)": {
                "formula": "C(8)H(6)N(6)O(2)"  # 219.0624999240
            },
            "cyclic N6-threonylcarbamoyladenosine -Ribose -C(4)H(7)N(1)O(2)": {
                "formula": "C(6)H(3)N(5)O(1)"  # 162.0410362034
            },
        },
    },
    "guanosine": {
        "references": [],
        "comments": "",
        "precursors": {
            "guanosine dimer": {"formula": "C(20)H(26)N(10)O(10)"},
            "guanosine": {"formula": "C(10)H(13)N(5)O(5)"},
            "guanosine K adduct": {"formula": "C(10)H(12)K(1)N(5)O(5)"},
        },
        "fragments": {
            # 'guanosine': {
            #   'formula': 'C(10)H(13)N(5)O(5)'
            # },
            "guanosine -Ribose": {"formula": "C(5)H(5)N(5)O(1)"},  # 152.0566862678
            "guanosine -Ribose -H(1)N(1) +O(1)": {
                "formula": "C(5)H(4)N(4)O(2)",  # 153.0407018516
                "hcd": True,
            },
            "guanosine -Ribose -H(3)N(1)": {
                "formula": "C(5)H(2)N(4)O(1)",  # 135.0301371672
                "hcd": True,
            },
        },
    },
    "cytidine": {
        "references": [],
        "comments": "",
        "precursors": {
            "cytidine dimer": {},
            "cytidine": {},
        },
        "fragments": {
            # 'cytidine': {
            #     'formula': 'C(9)H(13)N(3)O(5)'
            # },
            "cytidine -Ribose": {"formula": "C(4)H(5)N(3)O(1)"}  # 112.0505382598
        },
    },
    "3-methylcytidine": {
        # m3C
        "comments": "No dimer is formed! Methyl group blocks dimer formation?",
        "references": ["Jora et al. 2018"],
        "precursors": {
            #
            # '3-methylcytidine dimer' : { 'formula' : '' },
            # 'cytidine dimer semi-labeled' : { 'formula' : '' },
            "3-methylcytidine": {"formula": ""}
        },
        "fragments": {
            "3-methylcytidine -Ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            },
            "3-methylcytidine -Ribose -H(3)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "3-methylcytidine -Ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "3-methylcytidine -Methylated ribose": {
                "formula": "C(4)H(5)N(3)O(1)",  # 112.0505382598
                # 'hcd' : True
            },
            # '3-methylcytidine -Ribose' : {
            #     'formula': 'C(5)H(7)N(3)O(1)',
            #     'hcd': True
            # },
        },
        "dimer_possible": False,
    },
    "5-methylcytidine": {
        # m5C
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "5-methylcytidine dimer": {},
            "5-methylcytidine": {},
        },
        "fragments": {
            "5-methylcytidine -Ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            },
            "5-methylcytidine -Ribose -H(3)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "5-methylcytidine -Ribose -H(2)O(1)": {
                "formula": "C(5)H(5)N(3)",  # 108.0556236398
                "hcd": True,
            },
            # '5-methylcytidine -Ribose -C(1)H(1)N(1)O(1)' : {
            #     'formula': 'C(4)H(6)N(2)', # 83.0603746680
            #     'hcd': True
            # },
        },
        "exclusion_fragments": {
            "5-methylcytidine -Methylated ribose": {
                "formula": "C(4)H(5)N(3)O(1)",  # 112.0505382598
                "hcd": False,
            },
            "5-methylcytidine -Ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
        },
    },
    "N4-methylcytidine": {
        # m4C
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "N4-methylcytidine dimer": {},
            "N4-methylcytidine": {},
        },
        "fragments": {
            "N4-methylcytidine -Ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            },
            "N4-methylcytidine -Ribose -H(3)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "N4-methylcytidine -Ribose -H(2)O(1)": {
                "formula": "C(5)H(5)N(3)",  # 108.0556236398
                "hcd": True,
            },
            "N4-methylcytidine -Ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
            "N4-methylcytidine -Ribose -C(1)H(1)N(1)O(1)": {
                "formula": "C(4)H(6)N(2)",  # 83.0603746680
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "N4-methylcytidine -Methylated ribose": {
                "formula": "C(4)H(5)N(3)O(1)"  # 112.0505382598
            }
        },
    },
    "2′-O-methylcytidine": {
        # Cm
        "comments": "",
        "references": [],
        "precursors": {"2′-O-methylcytidine": {"formula": ""}},
        "fragments": {
            "2′-O-methylcytidine -Methylated ribose": {
                "formula": "C(4)H(5)N(3)O(1)"  # 112.0505382598
            },
            "2′-O-methylcytidine -Methylated ribose -H(3)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
            # '2′-O-methylcytidine -Methylated ribose -C(1)H(1)N(1)O(1)' : {
            #     'formula': 'C(3)H(4)N(2)', # 69.044
            #     'hcd': True
            # },
        },
        "exclusion_fragments": {
            "2′-O-methylcytidine -Ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            }
        },
    },
    "5-methoxycarbonylmethyluridine": {
        # mcm5U
        "comments": "",
        "references": [],
        "fragments": {
            "5-methoxycarbonylmethyluridine -Ribose": {
                "formula": "C(7)H(8)N(2)O(4)",  # 185.0556832124
                "hcd": False,
            },
            "5-methoxycarbonylmethyluridine -Ribose -C(1)H(4)O(1)": {
                "formula": "C(6)H(4)N(2)O(3)",  # 153.0294684636
                "hcd": False,
            },
            "5-methoxycarbonylmethyluridine -Ribose -C(2)H(4)O(2)": {
                "formula": "C(5)H(4)N(2)O(2)"  # 125.0345538436
            },
            "5-methoxycarbonylmethyluridine -Ribose -C(3)H(3)N(1)O(2)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            # '5-methoxycarbonylmethyluridine -Ribose -C(3)H(5)N(1)O(1)' : {
            #     'formula': 'C(4)H(3)N(1)O(1)', # 82.0287401874
            #     'hcd':True
            # },
            # '5-methoxycarbonylmethyluridine -Ribose -C(4)H(5)N(1)O(2)' : {
            #     'formula': 'C(3)H(3)N(1)', # 54.0338255674
            #     'hcd':True
            # },
        },
        "precursors": {
            "5-methoxycarbonylmethyluridine": {},
            "5-methoxycarbonylmethyluridine +H(3)N(1)": {
                "formula": "C(12)H(19)N(3)O(8)"
            },
            "5-methoxycarbonylmethyluridine dimer": {},
            "5-methoxycarbonylmethyluridine dimer +H(3)N(1)": {
                "formula": "C(24)H(35)N(5)O(16)"
            },
        },
    },
    "5-methoxycarbonylmethyl-2-thiouridine": {
        # mcm5s2U
        "comments": "",
        "references": [],
        "fragments": {
            "5-methoxycarbonylmethyl-2-thiouridine -Ribose": {
                "formula": "C(7)H(8)N(2)O(3)S(1)",  # 201.0328397664
                "hcd": False,
            },
            "5-methoxycarbonylmethyl-2-thiouridine -Ribose -C(1)H(4)O(1)": {
                "formula": "C(6)H(4)N(2)O(2)S(1)",  # 169.0066250176
                "hcd": False,
            },
            "5-methoxycarbonylmethyl-2-thiouridine -Ribose -C(2)H(4)O(2)": {
                "formula": "C(5)H(4)N(2)O(1)S(1)",  # 141.0117103976
                "hcd": True,
            },
            "5-methoxycarbonylmethyl-2-thiouridine -Ribose -C(3)H(3)N(1)O(1)S(1)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
        },
        "precursors": {"5-methoxycarbonylmethyl-2-thiouridine": {}},
    },
    "5-carbamoylmethyluridine": {
        # ncm5U
        "comments": "",
        "references": [],
        "precursors": {"5-carbamoylmethyluridine": {}},
        "fragments": {
            "5-carbamoylmethyluridine -Ribose": {
                "formula": "C(6)H(7)N(3)O(3)",  # 170.0560175642
                "hcd": False,
            },
            "5-carbamoylmethyluridine -Ribose -N(1)H(3)": {
                "formula": "C(6)H(4)N(2)O(3)",  # 153.0294684636
                "hcd": False,
            },
            "5-carbamoylmethyluridine -Ribose -C(1)H(3)N(1)O(1)": {
                "formula": "C(5)H(4)N(2)O(2)"  # 125.0345538436
            },
            "5-carbamoylmethyluridine -Ribose -C(2)H(2)N(2)O(1)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "5-carbamoylmethyluridine -Ribose -C(2)H(4)N(2)O(2)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
        },
    },
    "5-methoxycarbonylmethyl-2′-O-methyluridine": {
        # mcm5Um
        "comments": "",
        "references": [],
        "precursors": {
            "5-methoxycarbonylmethyl-2′-O-methyluridine": {
                "formula": "C(13)H(18)N(2)O(8)"  # 331.1135920144
            }
        },
        "fragments": {
            "5-methoxycarbonylmethyl-2′-O-methyluridine -Methylated ribose": {
                "formula": "C(7)H(8)N(2)O(4)",  # 185.0556832124
                "hcd": False,
            },
            "5-methoxycarbonylmethyl-2′-O-methyluridine -Methylated ribose -C(2)H(4)O(2)": {
                "formula": "C(5)H(4)N(2)O(2)"  # 125.0345538436
            },
            "5-methoxycarbonylmethyl-2′-O-methyluridine -Methylated ribose -C(3)H(3)N(1)O(2)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "5-methoxycarbonylmethyl-2′-O-methyluridine -Methylated ribose -C(3)H(4)N(1)O(3)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
        },
    },
    "5-carboxymethyluridine": {
        # cm5U
        "comments": "",
        "references": [],
        "precursors": {"5-carboxymethyluridine": {}},
        "fragments": {
            "5-carboxymethyluridine -Ribose": {
                "formula": "C(6)H(6)N(2)O(4)",  # 171.0400331480
                "hcd": False,
            },
            "5-carboxymethyluridine -Ribose -H(2)O(1)": {
                "formula": "C(6)H(4)N(2)O(3)",  # 153.0294684636
                "hcd": False,
            },  # 153
            "5-carboxymethyluridine -Ribose -C(1)H(2)O(2)": {
                "formula": "C(5)H(4)N(2)O(2)"  # 125.0345538436
            },
            "5-carbamoylmethyluridine -Ribose -C(2)H(2)N(2)O(1)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "5-carbamoylmethyluridine -Ribose -C(2)H(4)N(2)O(2)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
        },
    },
    "5-carbamoylmethyl-2-thiouridine": {
        # ncm5s2U
        "comments": "",
        "references": [],
        "precursors": {"5-carbamoylmethyl-2-thiouridine": {}},
        "fragments": {
            "5-carbamoylmethyl-2-thiouridine -Ribose": {
                "formula": "C(6)H(7)N(3)O(2)S(1)",  # 186.0331741182
                "hcd": False,
            },
            "5-carbamoylmethyl-2-thiouridine -Ribose -H(3)N(1)": {
                "formula": "C(6)H(4)N(2)O(2)S(1)",  # 169.0066250176
                "hcd": False,
            },
            "5-carbamoylmethyl-2-thiouridine -Ribose -C(1)H(3)O(1)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)S(1)"  # 141.0117103976
            },
            "5-carbamoylmethyl-2-thiouridine -Ribose -C(2)H(2)N(2)O(1)S(1)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "5-carbamoylmethyl-2-thiouridine -Ribose -C(2)H(4)N(2)O(2)S(1)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
        },
    },
    "5-carboxymethyl-2-thiouridine": {
        # cm5s2U
        "comments": "",
        "references": [],
        "precursors": {"5-carboxymethyl-2-thiouridine": {}},
        "fragments": {
            "5-carboxymethyl-2-thiouridine -Ribose": {
                "formula": "C(6)H(6)N(2)O(3)S(1)",  # 187.01718970197
                "hcd": False,
            },
            "5-carboxymethyl-2-thiouridine -Ribose -C(1)H(2)O(2)": {
                "formula": "C(5)H(4)N(2)O(1)S(1)"  # 141.0117103976
            },  # 141
            "5-carboxymethyl-2-thiouridine -Ribose -C(2)H(1)N(2)O(1)S(1)": {
                "formula": "C(4)H(5)N(1)O(2)",  # 100.0393048718
                "hcd": True,
            },
            "5-carboxymethyl-2-thiouridine -Ribose -C(2)H(3)N(1)O(2)S(1)": {
                "formula": "C(4)H(3)N(1)O(1)",  # 82.0287401874
                "hcd": True,
            },
        },
    },
    "5-carbamoylhydroxymethyluridine": {
        # nchm5U
        "comments": "",
        "references": [],
        "precursors": {"5-carbamoylhydroxymethyluridine": {}},
        "fragments": {
            "5-carbamoylhydroxymethyluridine -Ribose": {
                "formula": "C(6)H(7)N(3)O(4)",  # 186.0509321842
                "hcd": False,
            },
            "5-carbamoylhydroxymethyluridine -Ribose -C(1)H(3)N(1)O(1)": {
                "formula": "C(5)H(4)N(2)O(3)",  # 141.0294684636
                "hcd": True,
            },
            "5-carbamoylhydroxymethyluridine -Ribose -C(1)H(2)O(2)": {
                "formula": "C(5)H(5)N(3)O(2)",  # 140.0454528798
                "hcd": True,
            },
            "5-carbamoylhydroxymethyluridine -Ribose -C(1)H(2)O(2)": {
                "formula": "C(4)H(4)N(2)O(1)",  # 97.039639
                "hcd": True,
            },
            "5-carbamoylhydroxymethyluridine -Ribose -C(1)H(2)O(2)": {
                "formula": "C(3)H(3)N(1)O(1)",  # 70.0287401874
                "hcd": True,
            },
        },
    },
    "5-methyl-2-thiouridine": {
        # m5s2U
        "comments": "",
        "references": [],
        "precursors": {"5-methyl-2-thiouridine": {"formula": "C(10)H(14)N(2)O(5)S(1)"}},
        "fragments": {
            "5-methyl-2-thiouridine -Ribose": {
                "formula": "C(5)H(6)N(2)O(1)S(1)"  # 143.0273604620
            },
            "5-methyl-2-thiouridine -Ribose -N(1)H(3)": {
                "formula": "C(5)H(3)N(1)O(1)S(1)",  # 126.0008113614
                "hcd": True,
            },
        },
    },
    "5-aminomethyl-2-thiouridine": {
        # nm5s2U
        "comments": "",
        "references": [],
        "precursors": {
            "5-aminomethyl-2-thiouridine": {"formula": "C(10)H(15)N(3)O(5)S(1)"}
        },
        "fragments": {
            "5-aminomethyl-2-thiouridine -Ribose": {
                "formula": "C(5)H(7)N(3)O(1)S(1)"  # 158.038663833319
            },
            # "5-aminomethyl-2-thiouridine -Ribose -H(1)S(1)": {
            #     "formula": "C(5)H(6)N(3)O(1)",# 125.05836329197
            #     "hcd" : True
            #     not verified!!!
            # },
        },
    },
    "5-hydroxyuridine": {
        # ho5U
        "comments": "",
        "references": ["Nelson and McCloskey 1994"],
        "precursors": {"5-hydroxyuridine": {}},
        "fragments": {
            "5-hydroxyuridine -Ribose": {
                "formula": "C(4)H(4)N(2)O(3)"  # 129.0294684636
            },
            # "5-hydroxyuridine -Ribose -N(1)H(3)": {
            #     "formula": "C(4)H(1)N(1)O(3)", # 112.0029193630
            #     "hcd":True
            # },
            # "5-hydroxyuridine -Ribose -C(1)H(3)N(1)O(1)": {
            #     "formula": "C(3)H(1)N(1)O(2)", # 84.0080047430
            #     "hcd":True
            # },
        },
    },
    # "5,2′-O-dimethyluridine": {
    #     # m5Um
    #     "comments" : "Fragments similar to m5U (less similar to m3U, Verification required",
    #     "references" : ["Jora et al. 2018"],
    #     "precursors": {
    #         "5,2′-O-dimethyluridine": {
    #             "formula": "C(11)H(16)N(2)O(6)" # 273.1081127100
    #         }
    #     },
    #     "fragments": {
    #         "5,2′-O-dimethyluridine -Ribose": {
    #             "formula": "C(6)H(8)N(2)O(2)"  # 141.06585397237
    #         },
    #         "5,2′-O-dimethyluridine -Methylated ribose -N(1)H(3)": {
    #             "formula": "C(5)H(3)N(1)O(2)",  # 110.0236548074
    #             "hcd": True,
    #         },
    #         "5,2′-O-dimethyluridine -Methylated ribose -H(2)O(1)": {
    #             "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
    #             "hcd": True,
    #         },
    #         # '5,2′-O-dimethyluridine -Ribose -C(1)H(1)N(1)O(1)' : {
    #         #     'formula': 'C(4)H(5)N(1)O(1)', #84.0443902518
    #         #     'hcd':True
    #         # },
    #     },
    # },
    # "3,2′-O-dimethyluridine": {
    #     # m3Um
    #     "comments" : "Similarity to m3U assumed, Verfication required",
    #     "references" : ["Jora et al. 2018"],
    #     "precursors": {
    #         "3,2′-O-dimethyluridine": {
    #             "formula": "C(11)H(16)N(2)O(6)" # 273.1081127100
    #         }
    #     },
    #     "fragments": {
    #         "3,2′-O-dimethyluridine -Ribose": {
    #             "formula": "C(6)H(8)N(2)O(2)",  # 141.06585397237
    #             "hcd": True,
    #         },
    #         "3,2′-O-dimethyluridine -Ribose -C(1)H(5)N(1)": {
    #             "formula": "C(4)H(1)N(1)O(2)",  # 96.0080047430
    #             "hcd": True,
    #         },
    #     },
    # },
    "5-hydroxymethylcytidine": {
        "comments": "Elutes with m3C. Verification required. You et al. report peak 124.1",
        "references": ["You et al. 2019"],
        "precursors": {"5-hydroxymethylcytidine": {"formula": ""}},
        "fragments": {
            "5-hydroxymethylcytidine -Ribose": {
                "formula": "C(5)H(7)N(3)O(2)",  # 142.0611029442
                "hcd": False,
            },
            "5-hydroxymethylcytidine -Ribose -O(1)": {
                "formula": "C(5)H(7)N(3)O(1)",  #  126.0661883242
                "hcd": True,
            },
            "5-hydroxymethylcytidine -Ribose -H(3)N(1)O(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "5-hydroxymethylcytidine -Ribose -C(1)H(5)N(1)O(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
        },
    },
    "N4,2′-O-dimethylcytidine": {
        # m4Cm
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"N4,2′-O-dimethylcytidine": {}},
        "fragments": {
            "N4,2′-O-dimethylcytidine -Methylated ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            },
            "N4,2′-O-dimethylcytidine -Methylated ribose -H(3)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "N4,2′-O-dimethylcytidine -Methylated ribose -H(2)O(1)": {
                "formula": "C(5)H(5)N(3)",  # 108.0556236398
                "hcd": True,
            },
            "N4,2′-O-dimethylcytidine -Methylated ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
            "N4,2′-O-dimethylcytidine -Methylated ribose -C(1)H(1)N(1)O(1)": {
                "formula": "C(4)H(6)N(2)",  # 83.0603746680
                "hcd": True,
            },
        },
    },
    "5,2′-O-dimethylcytidine": {
        # m5Cm
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"5,2′-O-dimethylcytidine": {}},
        "fragments": {
            "5,2′-O-dimethylcytidine -Methylated ribose": {
                "formula": "C(5)H(7)N(3)O(1)"  # 126.0661883242
            },
            "5,2′-O-dimethylcytidine -Methylated ribose -H(3)N(1)": {
                "formula": "C(5)H(4)N(2)O(1)",  # 109.0396392236
                "hcd": True,
            },
            "5,2′-O-dimethylcytidine -Methylated ribose -H(2)O(1)": {
                "formula": "C(5)H(5)N(3)",  # 108.0556236398
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "5,2′-O-dimethylcytidine -Methylated ribose -C(1)H(5)N(1)": {  #
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            }
        },
    },
    "N4,N4-dimethylcytidine": {
        # m4Cm
        "comments": "Verification required",
        "references": [],
        "precursors": {"N4,N4-dimethylcytidine": {}},
        "fragments": {
            "N4,N4-dimethylcytidine -Ribose": {
                "formula": "C(6)H(9)N(3)O(1)"  # 140.0818383886
            }
        },
    },
    "1-methylguanosine": {
        # m1G
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"1-methylguanosine": {}, "1-methylguanosine -Ribose": {}},
        "fragments": {
            "1-methylguanosine -Ribose": {
                "formula": "C(6)H(7)N(5)O(1)"  # 166.0723363322
            },
            "1-methylguanosine -Ribose -H(1)N(1) +O(1)": {
                "formula": "C(6)H(6)N(4)O(2)",  # 167.0563519160
                "hcd": True,
            },
            "1-methylguanosine -Ribose -C(1)H(3)N(1) +O(1)": {
                "formula": "C(5)H(4)N(4)O(2)",  # 153.0407018516
                "hcd": True,
            },
            "1-methylguanosine -Ribose -H(3)N(1)": {
                "formula": "C(6)H(4)N(4)O(1)",  # 149.0457872316
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "7-O-methylguanosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)O(1)"  # 152.0566862678
            }
        },
    },
    "N2-methylguanosine": {
        # m2G
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"N2-methylguanosine": {}, "N2-methylguanosine -Ribose": {}},
        "fragments": {
            "N2-methylguanosine -Ribose": {
                "formula": "C(6)H(7)N(5)O(1)",  # 166.0723363322
                "hcd": False,
            },
            "N2-methylguanosine -Ribose -H(1)N(1) +O(1)": {
                "formula": "C(6)H(6)N(4)O(2)",  # 167.0563519160
                "hcd": True,
            },
            "N2-methylguanosine -Ribose -H(3)N(1)": {
                "formula": "C(6)H(4)N(4)O(1)",  # 149.0457872316
                "hcd": True,
            },
            "N2-methylguanosine -Ribose -C(2)H(2)N(2) +O(1)": {
                "formula": "C(4)H(5)N(3)O(2)",  # 128.0454528798
                "hcd": True,
            },
            "N2-methylguanosine -Ribose -C(2)H(4)N(2)": {
                "formula": "C(4)H(3)N(3)O(1)",  # 110.0348881954
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "N2-methylguanosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)O(1)"  # 152.0566862678
            },
            "N2-methylguanosine -Ribose": {
                "formula": "C(6)H(7)N(5)O(1)",  # 166.0723363322
                "hcd": True,
            },
        },
    },
    "7-methylguanosine": {
        # m7G
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"7-methylguanosine": {}, "7-methylguanosine -Ribose": {}},
        "fragments": {
            "7-methylguanosine -Ribose": {
                "formula": "C(6)H(7)N(5)O(1)"  # 166.0723363322
            },
            "7-methylguanosine -Ribose -H(1)N(1) +O(1)": {
                "formula": "C(6)H(6)N(4)O(2)",  # 167.0563519160
                "hcd": True,
            },
            "7-methylguanosine -Ribose -H(3)N(1)": {
                "formula": "C(6)H(4)N(4)O(1)",  # 149.0457872316
                "hcd": True,
            },
            # '7-methylguanosine -Ribose -C(1)N(2) +O(1)': { #too low abudant for 5% TIC threshold
            #     'formula': 'C(5)H(7)N(3)O(2)', # 142
            #     'hcd': True,
            # },
            "7-methylguanosine -Ribose -C(1)H(2)N(2)": {
                "formula": "C(5)H(5)N(3)O(1)",  # 124.0505382598
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "7-O-methylguanosine -Methylated ribose": {"formula": "C(5)H(5)N(5)O(1)"}
        },
    },
    "2′-O-methylguanosine": {
        # Cm
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "2′-O-methylguanosine": {},
            "2′-O-methylguanosine -Methylated ribose": {},
            "2′-O-methylguanosine K adduct": {},
            "2′-O-methylguanosine dimer": {},
            "2′-O-methylguanosine K adduct dimer": {},
        },
        "fragments": {
            "2′-O-methylguanosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)O(1)"  # 152.0566862678
            },
            "2′-O-methylguanosine -Methylated ribose -H(1)N(1) +O(1)": {
                "formula": "C(5)H(4)N(4)O(2)",  # 153.0407018516
                "hcd": True,
            },
            "2′-O-methylguanosine -Methylated ribose -H(3)N(1)": {
                "formula": "C(5)H(2)N(4)O(1)",  # 135.0301371672
                "hcd": True,
            },
            "2′-O-methylguanosine -Methylated ribose -C(1)H(2)N(2)": {
                "formula": "C(4)H(3)N(3)O(1)",  # 110.0348881954
                "hcd": True,
            },
            # '7-methylguanosine -Ribose -C(1)N(2) +O(1)': { #too low abudant for 5% TIC threshold
            #     'formula': 'C(5)H(7)N(3)O(2)', # 142
            #     'hcd': True,
            # },
        },
        "exclusion_fragments": {
            "2′-O-methylguanosine -Ribose": {
                "formula": "C(6)H(7)N(5)O(1)",  # 166.0723363322
            }
        },
    },
    "N2,N2-dimethylguanosine": {
        # m2,2G
        "comments": "",
        "references": ["Giessing et al. 2011"],
        "precursors": {
            "N2,N2-dimethylguanosine": {},
            "N2,N2-dimethylguanosine dimer": {},
        },
        "fragments": {
            "N2,N2-dimethylguanosine -Ribose": {
                "formula": "C(7)H(9)N(5)O(1)"  # 180.08798639657002
            },
            "N2,N2-dimethylguanosine -Ribose -C(2)H(5)N(1) +O(1)": {
                "formula": "C(5)H(4)N(4)O(2)",  # 153.0407018516
                "hcd": True,
            },
            "N2,N2-dimethylguanosine -Ribose -C(2)H(7)N(1)": {
                "formula": "C(5)H(2)N(4)O(1)",  # 135.0301371672
                "hcd": True,
            },
            "N2,N2-dimethylguanosine -Ribose -C(3)H(6)N(2)": {
                "formula": "C(4)H(3)N(3)O(1)",  # 110.0348881954
                "hcd": True,
            },
        },
    },
    # "7-aminomethyl-7-deazaguanosine": {
    #     # is not validated!
    #     "precursors": {"7-aminomethyl-7-deazaguanosine": {}},
    #     "fragments": {
    #         "7-aminomethyl-7-deazaguanosine -Ribose": {
    #             "formula": "C(7)H(9)N(5)O(1)"  # 180.0879863966
    #         },
    #         "7-aminomethyl-7-deazaguanosine -Ribose -H(3)N(1)": {
    #             "formula": "C(7)H(6)N(4)O(1)",  # 163.0614372960
    #             "hcd": True,
    #         },
    #     },
    # },
    "adenosine": {
        # A
        "comments": "",
        "references": [],
        "precursors": {"adenosine": {}},
        "fragments": {
            "adenosine -Ribose": {"formula": "C(5)H(5)N(5)"},  # 136.0617716478
            "adenosine -Ribose -H(3)N(1)": {
                "formula": "C(5)H(2)N(4)",  # 119.03522254717
                "hcd": True,
            },
        },
    },
    "2′-O-methyladenosine": {
        # Am
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "2′-O-methyladenosine": {"formula": "C(11)H(15)N(5)O(4)"},  # 282.1196804498
            "2′-O-methyladenosine dimer": {},
        },
        "fragments": {
            "2′-O-methyladenosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            },
            # '2′-O-methyladenosine': {'formula': 'C(11)H(15)N(5)O(4)'}
        },
        "exclusion_fragments": {
            "2′-O-methyladenosine -Ribose": {
                "formula": "C(6)H(7)N(5)"  # 150.0774217122
            }
        },
    },
    "1-methyladenosine": {
        # m1A
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "1-methyladenosine": {"formula": "C(11)H(15)N(5)O(4)"},
            "1-methyladenosine dimer": {},
            "1-methyladenosine K adduct dimer": {},
        },
        "fragments": {
            "1-methyladenosine -Ribose": {"formula": "C(6)H(7)N(5)"},  # 150.0774217122
            "1-methyladenosine -Ribose -N(1)H(3)": {
                "formula": "C(6)H(4)N(4)",  # 133.0508726116
                "hcd": True,
            },
            "1-methyladenosine -Ribose -C(2)N(1)H(3)": {
                "formula": "C(4)H(4)N(4)",  # 109.0508726116
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "1-methyladenosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            }
        },
    },
    "1,2′-O-dimethyladenosine": {
        # m1Am
        "comments": "",
        "references": [],
        "precursors": {"1,2′-O-dimethyladenosine": {"formula": "C(12)H(17)N(5)O(4)"}},
        "fragments": {
            "1,2′-O-dimethyladenosine -Methylated ribose": {
                "formula": "C(6)H(7)N(5)"  # 150.0774217122
            },
            "1,2′-O-dimethyladenosine -Methylated ribose -N(1)H(3)": {
                "formula": "C(6)H(4)N(4)",  # 133.0508726116
                "hcd": True,
            },
            "1,2′-O-dimethyladenosine -Methylated ribose -C(2)N(1)H(3)": {
                "formula": "C(4)H(4)N(4)",  # 109.0508726116
                "hcd": True,
            },
        },
    },
    "N6,N6-dimethyladenosine": {
        # m6,6A
        "comments": "",
        "references": ["Chan et al. 2011"],
        "precursors": {
            "N6,N6-dimethyladenosine": {
                "formula": "C(12)H(17)N(5)O(4)"  # 296.1353305142
            }
        },
        "fragments": {
            "N6,N6-dimethyladenosine -Ribose": {
                "formula": "C(7)H(9)N(5)"  # 164.0930717766
            },
            "N6,N6-dimethyladenosine -Ribose -C(1)H(3)": {
                "formula": "C(6)H(6)N(5)",  # 149.06959667997
                "hcd": True,
            },
            "N6,N6-dimethyladenosine -Ribose -C(2)H(5)N(1)": {
                "formula": "C(5)H(4)N(4)",  # 121.05087261157
                "hcd": True,
            },
        },
    },
    "N6,2′-O-dimethyladenosine": {
        # m6Am
        "comments": "",
        "references": [],
        "precursors": {
            "N6,2′-O-dimethyladenosine": {
                "formula": "C(12)H(17)N(5)O(4)"  # 296.1353305142
            }
        },
        "fragments": {
            "N6,2′-O-dimethyladenosine -Methylated ribose": {
                "formula": "C(6)H(7)N(5)"  # 150.0774217122
            },
            "N6,2′-O-dimethyladenosine -Methylated ribose -C(1)H(1)N(1)": {
                "formula": "C(5)H(6)N(4)",  # 123.0665226760
                "hcd": True,
            },
            "N6,2′-O-dimethyladenosine -Methylated ribose -C(2)H(4)N(2)": {
                "formula": "C(4)H(3)N(3)",  # 94.0399735754
                "hcd": True,
            },
        },
    },
    "2,8-dimethyladenosine": {
        # m2,8A
        "comments": "Verification required",
        "references": ["Giessing et al. 2009"],
        "precursors": {"2,8-dimethyladenosine": {}},
        "fragments": {
            "2,8-dimethyladenosine -Ribose": {
                "formula": "C(7)H(9)N(5)"  # 164.0930717766
            },
            "2,8-dimethyladenosine -Ribose -H(3)N(1)": {
                "formula": "C(7)H(6)N(4)",  # 147.0665226760
                "hcd": True,
            },
        },
    },
    "2-methyladenosine": {
        # m2A
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "2-methyladenosine": {},
            "2-methyladenosine dimer": {},
            "2-methyladenosine K adduct dimer": {},
        },
        "fragments": {
            "2-methyladenosine -Ribose": {"formula": "C(6)H(7)N(5)"},  # 150.0774217122
            "2-methyladenosine -Ribose -N(1)H(3)": {
                "formula": "C(6)H(4)N(4)",  # 133.0508726116
                "hcd": True,
            },
            "2-methyladenosine -Ribose -C(2)N(2)H(6)": {
                "formula": "C(4)H(1)N(3)",  # 92.0243235110
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "2-methyladenosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            }
        },
    },
    "8-methyladenosine": {
        # m8A
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "8-methyladenosine": {},
            "8-methyladenosine dimer": {},
            "8-methyladenosine K adduct dimer": {},
        },
        "fragments": {
            "8-methyladenosine -Ribose": {"formula": "C(6)H(7)N(5)"},  # 150.0774217122
            "8-methyladenosine -Ribose -N(1)H(3)": {
                "formula": "C(6)H(4)N(4)",  # 133.0508726116
                "hcd": True,
            },
            "8-methyladenosine -Ribose -C(1)N(2)H(2)": {
                "formula": "C(5)H(5)N(3)",  # 108.0556236398
                "hcd": True,
            },
            "8-methyladenosine -Ribose -C(1)N(2)H(4)": {
                "formula": "C(5)H(3)N(3)",  # 106.0399735754
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "8-methyladenosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            }
        },
    },
    "N6-methyladenosine": {
        # m6A
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "N6-methyladenosine": {},
            "N6-methyladenosine dimer": {},
            "N6-methyladenosine K adduct dimer": {},
        },
        "fragments": {
            "N6-methyladenosine -Ribose": {"formula": "C(6)H(7)N(5)"},  # 150.0774217122
            "N6-methyladenosine -Ribose -C(1)H(1)N(1)": {
                "formula": "C(5)H(6)N(4)",  # 123.0665226760
                "hcd": True,
            },
            "N6-methyladenosine -Ribose -C(2)H(4)N(2)": {
                "formula": "C(4)H(3)N(3)",  # 94.0399735754
                "hcd": True,
            },
        },
        "exclusion_fragments": {
            "N6-methyladenosine -Methylated ribose": {
                "formula": "C(5)H(5)N(5)"  # 136.0617716478
            }
        },
    },
    # 'N6-acetyladenosine': {
    #      "references":['Sauerwald et al. 2005'],
    #     'precursors':{
    #         'N6-acetyladenosine': {'formula': 'C(12)H(15)N(5)O(5)'}
    #     },
    #     'fragments' : {
    #         'N6-acetyladenosine -Ribose': {
    #             'formula': 'C(7)H(7)N(5)O(1)',
    #             'hcd' : False
    #         },
    #         'N6-acetyladenosine -Ribose -C(1)O(1)': {
    #             'formula': 'C(6)H(7)N(5)',
    #             'hcd' : True
    #         },
    #         'N6-acetyladenosine -Ribose -C(2)H(1)N(1)O(1)': {
    #             'formula': 'C(5)H(6)N(4)', #123.06
    #             'hcd' : True
    #         },
    #         'N6-acetyladenosine -Ribose -C(3)H(4)N(2)O(1)': {
    #             'formula': 'C(4)H(3)N(3)', #94.04
    #             'hcd' : True
    #         },
    #     },
    # },
    "N6-hydroxymethyladenosine": {
        # hm6A
        "comments": "Verification required",
        "references": ["modomics.org", "You et al. 2019", "Fu et al. 2013"],
        "fragments": {
            "N6-hydroxymethyladenosine": {
                "formula": "C(11)H(15)N(5)O(5)"  # 298.1145950698
            },
            "adenosine": {"formula": "C(10)H(13)N(5)O(4)"},  # 268.1040303854
            "adenosine -Ribose": {"formula": "C(5)H(5)N(5)"},  # 136.0617716478
        },
        "precursors": {"N6-hydroxymethyladenosine": {}},
    },
    "N4-acetylcytidine": {
        # ac4C
        "comments": "Verification required",
        "references": ["modomics.org"],
        "fragments": {
            "N4-acetylcytidine -Ribose": {
                "formula": "C(6)H(7)N(3)O(2)",  # 154.0611029442
                "hcd": False,
            },
            "cytidine -Ribose": {"formula": "C(4)H(5)N(3)O(1)"},  # 112.0505382598
            "cytidine -Ribose -C(1)H(5)N(1)": {
                "formula": "C(4)H(2)N(2)O(1)",  # 95.0239891592
                "hcd": True,
            },
        },
        "precursors": {"N4-acetylcytidine": {}, "N4-acetylcytidine K adduct": {}},
    },
    "5-formyl-2′-O-methylcytidine": {
        # f5Cm
        "comments": "Verification required",
        "references": ["modomics.org"],
        "fragments": {
            "5-formyl-2′-O-methylcytidine -Methylated ribose": {
                "formula": "C(5)H(5)N(3)O(2)"  # 140.0454528798
            }
        },
        "precursors": {"5-formyl-2′-O-methylcytidine": {}},
    },
    "4-thiouridine": {
        # s4U
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {"4-thiouridine": {}},
        "fragments": {
            "4-thiouridine -Ribose": {
                "formula": "C(4)H(4)N(2)O(1)S(1)"  # 129.0117103976
            },
            "4-thiouridine -Ribose -H(3)N(1)": {
                "formula": "C(4)H(1)N(1)O(1)S(1)",  #  111.9851612970
                "hcd": True,
            },
            "4-thiouridine -Ribose -C(1)H(1)N(1)O(1)": {
                "formula": "C(3)H(3)N(1)S(1)",  # 86.0058967414
                "hcd": True,
            },
        },
    },
    "2-thiouridine": {
        # s2U
        "comments": "",
        "references": ["Jora et al. 2018"],
        "precursors": {
            "2-thiouridine": {},
            "2-thiouridine isotope_1": {
                "formula": "+C(9)H(12)N(2)O(5)S(1)",
                "isotope_position": 1,
            },  # refers to pos 1 in most abundant precursors
            # default would be 0, following the ursgal nomenclature, see below
        },
        "fragments": {
            "2-thiouridine -Ribose": {
                "formula": "C(4)H(4)N(2)O(1)S(1)"  # 129.0117103976
            },
            "2-thiouridine -Ribose -H(3)N(1)": {
                "formula": "C(4)H(1)N(1)O(1)S(1)",  #  111.9851612970
                "hcd": True,
            },
            "2-thiouridine -Ribose -C(1)H(1)N(1)S(1)": {
                "formula": "C(3)H(3)N(1)O(1)",  #  70.0287401874
                "hcd": True,
            },
        },
    },
    "spiroiminodihydantoin": {
        "comments": "",
        "references": ["Martinez et al. 2007"],
        "precursors": {
            "spiroiminodihydantoin": {"formula": "C(10)H(13)N(5)O(7)"}  # 316.0887742454
        },
        "fragments": {
            "spiroiminodihydantoin -Ribose": {
                "formula": "C(5)H(5)N(5)O(3)",  # 184.0465155078
                "hcd": False,
            },
            "spiroiminodihydantoin -Ribose -C(2)H(2)N(2)O(1)": {
                "formula": "C(3)H(3)N(3)O(2)",  # 114.0298028154
                "hcd": True,
            },
            "spiroiminodihydantoin -Ribose -C(2)H(1)N(1)O(2)": {
                "formula": "C(3)H(4)N(4)O(1)",  # 113.0457872316
                "hcd": True,
            },
            "spiroiminodihydantoin -Ribose -C(3)H(2)N(2)O(2)": {
                "formula": "C(2)H(3)N(3)O(1)",  # 86.0348881954
                "hcd": True,
            },
        },
    },
    "queuosine": {
        # Q
        "comments": "",
        "references": ["Zallot et al. 2014", "Giessing et al. 2009"],
        "precursors": {"queuosine": {"formula": "C(17)H(23)N(5)O(7)"}},
        "fragments": {
            "queuosine -C(5)H(9)N(1)O(2)": {
                "formula": "C(12)H(14)N(4)O(5)",  # 295.10369603357
                "hcd": False,
            },
            "queuosine -C(10)H(17)N(1)O(6)": {
                "formula": "C(7)H(6)N(4)O(1)",  # 163.0614372959700
                # 'hcd' : False
            },
            "queuosine -C(10)H(20)N(2)O(6)": {
                "formula": "C(7)H(3)N(3)O(1)",  # 146.0348881954
                "hcd": True,
            },
            "queuosine -C(11)H(19)N(3)O(6)": {
                "formula": "C(6)H(4)N(2)O(1)",  # 121.0396392236
                "hcd": True,
            },
        },
    },
}
