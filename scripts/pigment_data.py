#!/usr/bin/env python3

# note : usage is pigment_data[ligand_code][state][variable]
pigment_data = {
        'CLA': { # note that chlorophyll Qy := S1, Qx := S2
            'D' : 4.0,
            'S1': {
                'energy': 14900.0,
                'lifetime': 4.0,
                'reorg': 0.0, # actually will need to read this in
                },
            'S2': {
                # 'energy': 16300.0,
                'energy': 15900.0,
                # 'lifetime': 4.0, # probably untrue
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'CHL': { # chlorophyll b
            'D' : 3.4,
            'S1': {
                'energy': 15300.0,
                'lifetime': 4.0,
                'reorg': 0.0, # actually will need to read this in
                },
            'S2': {
                'energy': 16700.0,
                # 'lifetime': 4.0, # probably untrue
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'KC1': {
            'S1': {
                # 'energy': 15970.0,
                'energy': 15600.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                },
            'S2': {
                # 'energy': 17330.0,
                'energy': 16700.0,
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'KC2': {
            'S1': {
                'energy': 15600.0,
                'lifetime': 4.0,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 16700.0,
                'lifetime': 0.0001, # probably untrue
                'reorg': 0.0,
                }
            },
        'A86': {
            'S1': {
                'energy': 15870.0,
                'lifetime': 0.01,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006, # 60fs
                'reorg': 0.0,
                }
            },
        'DD6': {
            'S1': {
                'energy': 14900.0,
                'lifetime': 0.01,
                'reorg': 0.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'LUT': {
            'D' : 0.1, # kieran carotenoid pathway table 1
            'S1': {
                'energy': 15400.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'NEX': {
            'D' : 0.15,
            'S1': { # no real data for this or XAT
                'energy': 15200.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        'XAT': {
            'D' : 0.18,
            'S1': { # no real data for this or XAT
                'energy': 15200.0,
                'lifetime': 0.01,
                'reorg': 450.0,
                },
            'S2': {
                'energy': 20000.0,
                'lifetime': 0.00006,
                'reorg': 0.0,
                }
            },
        }

# check where these are from!!!
# muh and renger
site_energies = {
        'CHL601' : 15405,
        'CLA602' : 14940,
        'CLA603' : 14850,
        'CLA604' : 14820,
        'CHL605' : 15465,
        'CHL606' : 15385,
        'CHL607' : 15225,
        'CHL608' : 15215,
        'CHL609' : 15475,
        'CLA610' : 14790,
        'CLA611' : 14950,
        'CLA612' : 14940,
        'CLA613' : 14840,
        'CLA614' : 14940,
        'LUT620' : 14525, # testing
        'LUT621' : 14475,
        # 'LUT620' : 14500, # guess from chris - both luteins 14500
        # 'LUT621' : 14500,
        # 'LUT620' : 19538, # 620/621 from kieran carotenoid paper
        # 'LUT621' : 19748,
        }

novod_energies = {
        'CHL601' : 15890,
        'CLA602' : 15160,
        'CLA603' : 15283,
        'CLA604' : 15460,
        'CHL605' : 15679,
        'CHL606' : 15851,
        'CHL607' : 15712,
        'CHL608' : 15763,
        'CHL609' : 15721,
        'CLA610' : 15073,
        'CLA611' : 15115,
        'CLA612' : 15097,
        'CLA613' : 15175,
        'CLA614' : 15264,
        }



vangelis_initial = {
        'CHL1601': 'CHL262',
        'CHL1605': 'CHL264',
        'CHL1606': 'CHL264',
        'CHL1607': 'CHL265',
        'CHL1608': 'CHL264',
        'CHL1609': 'CHL265',
        'CHL5601': 'CHL270',
        'CHL5605': 'CHL260',
        'CHL5606': 'CHL261',
        'CHL5607': 'CHL271',
        'CHL5608': 'CHL266',
        'CHL5609': 'CHL273',
        'CHL6601': 'CHL266',
        'CHL6605': 'CHL258',
        'CHL6606': 'CHL259',
        'CHL6607': 'CHL267',
        'CHL6608': 'CHL258',
        'CHL6609': 'CHL269',
        'CLA1602': 'CLA232',
        'CLA1603': 'CLA233',
        'CLA1604': 'CLA250',
        'CLA1610': 'CLA234',
        'CLA1611': 'CLA236',
        'CLA1612': 'CLA236',
        'CLA1613': 'CLA237',
        'CLA1614': 'CLA251',
        'CLA5602': 'CLA244',
        'CLA5603': 'CLA245',
        'CLA5604': 'CLA254',
        'CLA5610': 'CLA246',
        'CLA5611': 'CLA246',
        'CLA5612': 'CLA246',
        'CLA5613': 'CLA249',
        'CLA5614': 'CLA247',
        'CLA6602': 'CLA239',
        'CLA6603': 'CLA239',
        'CLA6604': 'CLA237',
        'CLA6610': 'CLA252',
        'CLA6611': 'CLA252',
        'CLA6612': 'CLA252',
        'CLA6613': 'CLA243',
        'CLA6614': 'CLA253',
        'LUT1620': 'LUT274',
        'LUT1621': 'LUT275',
        'LUT5620': 'LUT278',
        'LUT5621': 'LUT279',
        'LUT6620': 'LUT276',
        'LUT6621': 'LUT277',
        'NEX1623': 'NEX280',
        'NEX5623': 'NEX282',
        'NEX6623': 'NEX281',
        'XAT1622': 'XAT284',
        'XAT5622': 'XAT283',
        'XAT6622': 'XAT285',
        }

# if the keys are vangelis' residues then the values will be the crystal structure - in that case value[3] tells us which monomer we're in
# we get the key from "{}{}".format(arr[0][3], str(arr[0][4]))
