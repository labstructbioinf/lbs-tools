''' functions for parsing PLIP output interaction data '''
''' data parsed with this functions can be stored as pickle data '''

def parse_hbond(hb):
    ''' parse hydrogen bond type interaction, example of hbond class data from PLIP below '''
    ''' hbond(a=<pybel.Atom object at 0x7f2bbc16da90>, a_orig_idx=464, d=<pybel.Atom object at 0x7f2bbbcb1d10>,
            d_orig_idx=4916, h=<pybel.Atom object at 0x7f2bbcdaf1d0>, distance_ah=2.0442823843932194,
            distance_ad=2.5851673059978126, angle=113.21198137185546, type='strong', protisdon=False,
            resnr=61, restype='SER', reschain='A', resnr_l=400, restype_l='NAD',
            reschain_l='A', sidechain=True, atype='O3', dtype='O3') '''

    a = hb.a.type + ', ' + str(hb.a)
    d = hb.d.type + ', ' + str(hb.d)
    h = hb.h.type + ', ' + str(hb.h)

    record = {'a_orig_idx':hb.a_orig_idx, 'd_orig_idx':hb.d_orig_idx, 'distance_ah':hb.distance_ah,
              'distance_ad':hb.distance_ad, 'angle':hb.angle, 'type':hb.type, 'protisdon':hb.protisdon,
              'resnr':hb.resnr, 'restype':hb.restype, 'reschain':hb.reschain, 'resnr_l':hb.resnr_l,
              'restype_l':hb.restype_l, 'reschain_l':hb.reschain_l, 'sidechain':hb.sidechain,
              'atype':hb.atype, 'dtype':hb.dtype, 'a':a, 'd':d, 'h':h, 'inter_type':'hbond'}

    return record

def parse_hydroph_interaction(hi):
    ''' parse hydrophobic type interaction, example of hydroph_interaction class data from PLIP below '''
    ''' hydroph_interaction(bsatom=<pybel.Atom object at 0x7f2bbc146ed0>, bsatom_orig_idx=945,
            ligatom=<pybel.Atom object at 0x7f2bbbcb7d50>, ligatom_orig_idx=4948,
            distance=3.89281106656873, restype='LEU', resnr=122, reschain='A', restype_l='NAD',
            resnr_l=400, reschain_l='A') '''

    bsatom  = hi.bsatom.type + ', ' + str(hi.bsatom)
    ligatom = hi.ligatom.type + ', ' + str(hi.ligatom)

    record = {'bsatom_orig_idx':hi.bsatom_orig_idx, 'ligatom_orig_idx':hi.ligatom_orig_idx,
              'distance':hi.distance, 'restype':hi.restype, 'resnr':hi.resnr, 'reschain':hi.reschain,
              'restype_l':hi.restype_l, 'resnr_l':hi.resnr_l, 'reschain_l':hi.reschain_l, 'bsatom':bsatom,
              'ligatom':ligatom, 'inter_type':'hydroph_interaction'}

    return record

def parse_waterbridge(wb):
    ''' parse waterbridge type interaction, example of waterbridge class data from PLIP below '''
    ''' waterbridge(a=<pybel.Atom object at 0x7f2bbc163510>, a_orig_idx=1083, atype='O2',
            d=<pybel.Atom object at 0x7f2bbbcb78d0>, d_orig_idx=4939, dtype='O3',
            h=<pybel.Atom object at 0x7f2bbcdaf6d0>, water=<pybel.Atom object at 0x7f2bbbcd1190>,
            water_orig_idx=5100, distance_aw=3.9570770778442044, distance_dw=2.7366486804118657,
            d_angle=135.58856241979237, w_angle=112.18800830703512, type='first_deg', resnr=142,
            restype='ASN', reschain='A', resnr_l=400, restype_l='NAD', reschain_l='A', protisdon=False) '''

    a     = wb.a.type + ', ' + str(wb.a)
    d     = wb.d.type + ', ' + str(wb.d)
    h     = wb.h.type + ', ' + str(wb.h)
    water = wb.water.type + ', ' + str(wb.water)

    record = {'a_orig_idx':wb.a_orig_idx, 'atype':wb.atype, 'd_orig_idx':wb.d_orig_idx,
              'dtype':wb.dtype, 'water_orig_idx':wb.water_orig_idx, 'distance_aw':wb.distance_aw,
              'distance_dw':wb.distance_dw, 'd_angle':wb.d_angle, 'w_angle':wb.w_angle, 'type':wb.type,
              'resnr':wb.resnr, 'restype':wb.restype, 'reschain':wb.reschain, 'resnr_l':wb.resnr_l,
              'restype_l':wb.restype_l, 'reschain_l':wb.reschain_l, 'protisdon':wb.protisdon, 'a':a,
              'd':d, 'h':h, 'water':water, 'inter_type':'waterbridge'}

    return record

def parse_saltbridge(sb):
    ''' parse saltbridge type interaction, example of saltbridge class data from PLIP below '''
    ''' saltbridge(positive=pcharge(atoms=[<pybel.Atom object at 0x7f3d316b8650>, <pybel.Atom object at 0x7f3d316b86d0>,
                                           <pybel.Atom object at 0x7f3d316b8790>], atoms_orig_idx=[118, 120, 121],
                                           type='positive', center=[43.640666666666675, 13.101333333333335, -24.508666666666667],
                                           restype='ARG', resnr=23, reschain='A'),
                   negative=lcharge(atoms=[<pybel.Atom object at 0x7f3d311ea9d0>, <pybel.Atom object at 0x7f3d311f94d0>,
                                           <pybel.Atom object at 0x7f3d311f9a10>, <pybel.Atom object at 0x7f3d311f9c10>,
                                           <pybel.Atom object at 0x7f3d311f9e50>],
                                           orig_atoms=[<pybel.Atom object at 0x7f3d311f9e90>, <pybel.Atom object at 0x7f3d311f9f90>,
                                           <pybel.Atom object at 0x7f3d311f9c50>, <pybel.Atom object at 0x7f3d31216450>,
                                           <pybel.Atom object at 0x7f3d31216310>, <pybel.Atom object at 0x7f3d31216350>],
                                           atoms_orig_idx=[2599, 2599, 2610, 2611, 2604, 2607], type='negative',
                                           center=(41.562, 14.363, -21.049), fgroup='phosphate'),
                   distance=4.2287056333272215, protispos=True, resnr=23, restype='ARG', reschain='A', resnr_l=334,
                   restype_l='TLO', reschain_l='A') '''

    positive_type = str(type(sb.positive)).strip('<>\'').split('.')[-1]
    negative_type = str(type(sb.negative)).strip('<>\'').split('.')[-1]

    if positive_type == 'pcharge':
        positive = {'atoms':[ a.type + ', ' + str(a) for a in sb.positive.atoms ], 'type':sb.positive.type,
                    'center':sb.positive.center, 'restype':sb.positive.restype, 'resnr':sb.positive.resnr,
                    'reschain':sb.positive.reschain}
    elif positive_type == 'lcharge':
        positive = {'atoms':[ a.type + ', ' + str(a) for a in sb.positive.atoms ],
                    'orig_atoms':[ a.type + ', ' + str(a) for a in sb.positive.orig_atoms ],
                    'atoms_orig_idx':sb.positive.atoms_orig_idx, 'type':sb.positive.type,
                    'center':sb.positive.center, 'fgroup':sb.positive.fgroup}

    if negative_type == 'lcharge':
        negative = {'atoms':[ a.type + ', ' + str(a) for a in sb.negative.atoms ],
                    'orig_atoms':[ a.type + ', ' + str(a) for a in sb.negative.orig_atoms ],
                    'atoms_orig_idx':sb.negative.atoms_orig_idx, 'type':sb.negative.type,
                    'center':sb.negative.center, 'fgroup':sb.negative.fgroup}
    elif negative_type == 'pcharge':
        negative = {'atoms':[ a.type + ', ' + str(a) for a in sb.negative.atoms ], 'type':sb.negative.type,
                    'center':sb.negative.center, 'restype':sb.negative.restype, 'resnr':sb.negative.resnr,
                    'reschain':sb.negative.reschain}

    record = {'distance':sb.distance, 'protispos':sb.protispos, 'resnr':sb.resnr, 'restype':sb.restype,
              'reschain':sb.reschain, 'resnr_l':sb.resnr_l, 'restype_l':sb.restype_l, 'reschain_l':sb.reschain_l,
              'positive':positive, 'negative':negative, 'inter_type':'saltbridge'}

    return record

def parse_metal_complex(mc):
    ''' parse metal complex type interaction, example of metal_complex class data from PLIP below '''
    ''' metal_complex(metal=<pybel.Atom object at 0x7f3984d7c7d0>, metal_orig_idx=5430, metal_type='Fe',
                             target=metal_binding(atom=<pybel.Atom object at 0x7f3984b0b950>, atom_orig_idx=4742,
                             type='N', restype='HIS', resnr=92, reschain='D', location='protein.sidechain'),
                      target_orig_idx=4742, target_type='N', coordination_num=5, distance=2.148810135865895,
                      resnr=92, restype='HIS', reschain='D', restype_l='HEM', reschain_l='D', resnr_l=1147,
                      location='protein.sidechain', rms=20.264282843925933, geometry='square.pyramidal',
                      num_partners=5, complexnum=1) '''

    metal = mc.metal.type + ', ' + str(mc.metal)
    atom  = mc.target.atom.type + ', ' + str(mc.target.atom)

    target = {'atom':atom, 'atom_orig_idx':mc.target.atom_orig_idx, 'type':mc.target.type, 'restype':mc.target.restype,
              'resnr':mc.target.resnr, 'reschain':mc.target.reschain,'location':mc.target.location}

    record = {'metal_orig_idx':mc.metal_orig_idx, 'metal_type':mc.metal_type, 'target':target,
             'target_orig_idx':mc.target_orig_idx, 'target_type':mc.target_type,
              'coordination_num':mc.coordination_num,'distance':mc.distance, 'resnr':mc.resnr,
              'restype':mc.restype, 'reschain':mc.reschain, 'restype_l':mc.restype_l, 'reschain_l':mc.reschain_l,
              'resnr_l':mc.resnr_l, 'location':mc.location, 'rms':mc.rms, 'geometry':mc.geometry,
              'num_partners':mc.num_partners, 'complexnum':mc.complexnum, 'inter_type':'metal_complex'}

    return record

def parse_pistack(ps):
    ''' parse pi-stack type interaction, example of pistack class data from PLIP below '''
    ''' pistack(proteinring=aromatic_ring(atoms=[<pybel.Atom object at 0x7f3ce3ca4a90>, <pybel.Atom object at 0x7f3ce3ca4b10>,
                                                 <pybel.Atom object at 0x7f3ce3ca4b90>, <pybel.Atom object at 0x7f3ce3ca4c10>,
                                                 <pybel.Atom object at 0x7f3ce3ca4c90>, <pybel.Atom object at 0x7f3ce3ca4d10>],
                                          orig_atoms=[<pybel.Atom object at 0x7f3ce339f2d0>, <pybel.Atom object at 0x7f3ce339f050>,
                                                      <pybel.Atom object at 0x7f3ce339f3d0>, <pybel.Atom object at 0x7f3ce339f150>,
                                                      <pybel.Atom object at 0x7f3ce339f4d0>, <pybel.Atom object at 0x7f3ce339f250>],
                                          atoms_orig_idx=[189, 190, 191, 192, 193, 194], normal=array([-0.25022486,  0.47387021, -0.8442953 ]),
                                          obj=<openbabel.OBRing; proxy of <Swig Object of type 'OpenBabel::OBRing *' at 0x7f3ce33a3240> >,
                                          center=[-11.847666666666669, 18.087166666666665, -20.311999999999998], type='6-membered'),
                            ligandring=aromatic_ring(atoms=[<pybel.Atom object at 0x7f3ce3426f10>, <pybel.Atom object at 0x7f3ce3426f90>,
                                                            <pybel.Atom object at 0x7f3ce342c050>, <pybel.Atom object at 0x7f3ce342c250>,
                                                            <pybel.Atom object at 0x7f3ce342c2d0>, <pybel.Atom object at 0x7f3ce342c350>],
                                                     orig_atoms=[<pybel.Atom object at 0x7f3ce34319d0>, <pybel.Atom object at 0x7f3ce3431b10>,
                                                                 <pybel.Atom object at 0x7f3ce3431b90>, <pybel.Atom object at 0x7f3ce3431c10>,
                                                                 <pybel.Atom object at 0x7f3ce3431c90>, <pybel.Atom object at 0x7f3ce3431d10>],
                                                     atoms_orig_idx=[4558, 4559, 4560, 4564, 4565, 4566], normal=array([-0.26086629,  0.40432458, -0.87662444]),
                                                     obj=<openbabel.OBRing; proxy of <Swig Object of type 'OpenBabel::OBRing *' at 0x7f3ce3424cf0> >,
                                                     center=[-12.282166666666669, 16.281333333333333, -17.16933333333333], type='6-membered'),
                            distance=3.6505038084565196, angle=4.437374472274115, offset=1.3280892590992437, type='P', restype='TYR', resnr=22, reschain='A',
                            restype_l='NAD', resnr_l=701, reschain_l='A') '''

    proteinring = {'atoms':[ a.type + ', ' + str(a) for a in ps.proteinring.atoms ],
                   'orig_atoms':[ a.type + ', ' + str(a) for a in ps.proteinring.orig_atoms ],
                   'atoms_orig_idx':ps.proteinring.atoms_orig_idx, 'normal':ps.proteinring.normal,
                   'center':ps.proteinring.center, 'type':ps.proteinring.type}

    ligandring = {'atoms':[ a.type + ', ' + str(a) for a in ps.ligandring.atoms ],
                   'orig_atoms':[ a.type + ', ' + str(a) for a in ps.ligandring.orig_atoms ],
                   'atoms_orig_idx':ps.ligandring.atoms_orig_idx, 'normal':ps.ligandring.normal,
                   'center':ps.ligandring.center, 'type':ps.ligandring.type}

    record = {'distance':ps.distance, 'angle':ps.angle,'offset':ps.offset, 'type':ps.type, 'restype':ps.restype,
              'resnr':ps.resnr, 'reschain':ps.reschain, 'restype_l':ps.restype_l,'resnr_l':ps.resnr_l,
              'reschain_l':ps.reschain_l, 'proteinring':proteinring, 'ligandring':ligandring, 'inter_type':'pistack'}

    return record

def parse_pication(pc):
    ''' parse pi-cation type interaction, example of pication class data from PLIP below '''
    ''' pication(ring=aromatic_ring(atoms=[<pybel.Atom object at 0x7f3a4498ee50>, <pybel.Atom object at 0x7f3a4498eed0>,
                                           <pybel.Atom object at 0x7f3a4498ef50>, <pybel.Atom object at 0x7f3a4498efd0>,
                                           <pybel.Atom object at 0x7f3a4497f090>],
                                    orig_atoms=[<pybel.Atom object at 0x7f3a4497fd50>, <pybel.Atom object at 0x7f3a4497fdd0>,
                                                <pybel.Atom object at 0x7f3a4497fe50>, <pybel.Atom object at 0x7f3a4497fed0>,
                                                <pybel.Atom object at 0x7f3a4497ff50>],
                                    atoms_orig_idx=[56019, 56020, 56021, 56022, 56023], normal=array([-0.51123067, -0.20248671,  0.83524986]),
                                    obj=<openbabel.OBRing; proxy of <Swig Object of type 'OpenBabel::OBRing *' at 0x7f3a449f25d0> >,
                                    center=[71.6358, -36.4092, 173.539], type='5-membered'),
                 charge=pcharge(atoms=[<pybel.Atom object at 0x7f3a4466ad10>, <pybel.Atom object at 0x7f3a4466ad50>,
                                       <pybel.Atom object at 0x7f3a4466ae10>], atoms_orig_idx=[15283, 15285, 15286],
                                       type='positive', center=[69.00333333333334, -35.791, 176.72433333333333],
                                       restype='ARG', resnr=20, reschain='G'),
                 distance=4.178325075380756, offset=1.5475537726778732, type='regular', restype='ARG', resnr=20,
                 reschain='G', restype_l='CLA', resnr_l=104, reschain_l='G', protcharged=True) '''

    charge_type = str(type(pc.charge)).strip('<>\'').split('.')[-1]

    ring = {'atoms':[ a.type + ', ' + str(a) for a in pc.ring.atoms ],
            'orig_atoms':[ a.type + ', ' + str(a) for a in pc.ring.orig_atoms ],
            'atoms_orig_idx':pc.ring.atoms_orig_idx, 'normal':pc.ring.normal,
            'center':pc.ring.center, 'type':pc.ring.type}

    if charge_type == 'pcharge':
        charge = {'atoms':[ a.type + ', ' + str(a) for a in pc.charge.atoms ], 'type':pc.charge.type,
                  'center':pc.charge.center, 'restype':pc.charge.restype, 'resnr':pc.charge.resnr,
                  'reschain':pc.charge.reschain}
    elif charge_type == 'lcharge':
         charge = {'atoms':[ a.type + ', ' + str(a) for a in pc.charge.atoms ], 'orig_atoms':[ a.type + ', ' + str(a) for a in pc.charge.orig_atoms ],
                   'atoms_orig_idx':pc.charge.atoms_orig_idx, 'type':pc.charge.type, 'center':pc.charge.center, 'fgroup':pc.charge.fgroup}

    record = {'distance':pc.distance, 'offset':pc.offset, 'type':pc.type, 'restype':pc.restype, 'resnr':pc.resnr,
              'reschain':pc.reschain, 'restype_l':pc.restype_l, 'resnr_l':pc.resnr_l, 'reschain_l':pc.reschain_l,
              'protcharged':pc.protcharged, 'ring':ring, 'charge':charge, 'inter_type':'pication'}

    return record

def parse_halogenbond(ha):
    ''' parse halogen bond type interaction, example of halogenbond class data from PLIP below '''
    ''' halogenbond(acc=hal_acceptor(o=<pybel.Atom object at 0x7f46bf453a10>, o_orig_idx=711, y=<pybel.Atom object at 0x7f46bdf033d0>, y_orig_idx=710),
    acc_orig_idx=711, don=hal_donor(x=<pybel.Atom object at 0x7f46bc429390>, orig_x=<pybel.Atom object at 0x7f46bef9e550>, x_orig_idx=2468,
    c=<pybel.Atom object at 0x7f46bef9e9d0>, c_orig_idx=[2449]),
    don_orig_idx=2468, distance=3.1128318297010433, don_angle=177.25365384843246, acc_angle=127.80617506353835, restype='VAL', resnr=127,
    reschain='A', restype_l='3BM', resnr_l=1, reschain_l='A', donortype='I', acctype='O2', sidechain=False) '''

    o = ha.acc.o.type + ', ' + str(ha.acc.o)
    y = ha.acc.y.type + ', ' + str(ha.acc.y)
    x = ha.don.x.type + ', ' + str(ha.don.x)
    orig_x = ha.don.orig_x.type + ', ' + str(ha.don.orig_x)
    c = ha.don.c.type + ', ' + str(ha.don.c)

    hal_acceptor = {'o':o, 'o_orig_idx':ha.acc.o_orig_idx, 'y':y, 'y_orig_idx':ha.acc.y_orig_idx}
    hal_donor    = {'x':x, 'orig_x':orig_x, 'x_orig_idx':ha.don.x_orig_idx, 'c':c, 'c_orig_idx':ha.don.c_orig_idx}


    record = {'acc':hal_acceptor, 'acc_orig_idx':ha.acc_orig_idx, 'don':hal_donor, 'don_orig_idx':ha.don_orig_idx, 'distance':ha.distance, 'don_angle':ha.don_angle,
      'acc_angle':ha.acc_angle, 'restype':ha.restype, 'resnr':ha.resnr, 'reschain':ha.reschain, 'restype_l':ha.restype_l, 'resnr_l':ha.resnr_l,
      'reschain_l':ha.reschain_l, 'donortype':ha.donortype, 'acctype':ha.acctype, 'sidechain':ha.sidechain, 'inter_type':'halogenbond'}

    return record

def parse_inter(inter):
    inter_type = str(type(inter)).strip('<>\'').split('.')[-1]
    if inter_type == 'saltbridge':
        return parse_saltbridge(inter)
    elif inter_type == 'hydroph_interaction':
        return parse_hydroph_interaction(inter)
    elif inter_type == 'waterbridge':
        return parse_waterbridge(inter)
    elif inter_type == 'hbond':
        return parse_hbond(inter)
    elif inter_type == 'metal_complex':
        return parse_metal_complex(inter)
    elif inter_type == 'pistack':
        return parse_pistack(inter)
    elif inter_type == 'pication':
        return parse_pication(inter)
    elif inter_type == 'halogenbond':
        return parse_halogenbond(inter)
    else:
        print('WARNING: unrecognized interaction ', inter_type)
        return None
