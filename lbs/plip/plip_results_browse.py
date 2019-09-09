import pickle

plip_data = pickle.load(open('/home/kszczepaniak/Data/PiPred/Figure/PLIP_Figure/plip_interactions_all.pickle.txt', 'rb'))

train_pdb_scan_6 = pickle.load(open("/home/kszczepaniak/Data/PiPred/Figure/PLIP_Figure/id_data_clean_dict_FINAL_FINAL_6.p2"))

pdb_id = '4dja_A'

for p in plip_data.items():
    for i in p[1]:
        if i['restype_l'] == 'FAD':
            #print i
            print p[0], i['inter_type']
            if p[0] == pdb_id:
                print i

# HEM - od Staszka
# CLA - chlorophyl

# 4xk8_G: CLA, hp+pi-cat
# 2hun_B: NAD, hbond
# 2eq6_B: FAD, hbond
# 3slk_B: NDP, hbond (NADPH)


print train_pdb_scan_6[pdb_id][0]
print train_pdb_scan_6[pdb_id][1]
