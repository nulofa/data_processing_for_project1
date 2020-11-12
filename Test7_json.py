import numpy as np
import json

sam_cond = np.load('./sample_cond_dels2.npy')
sam_d = np.load('./sam_details_dels2.npy', allow_pickle=True)
sam_d=sam_d.item()
data_2d = np.load('./data2d_dels2.npy')
data_json = []
i=0
for sc in sam_cond:
    gse = sc.split('_')[0]

    sim = ''

    if sc in sam_d:
        Accesion_number,Title, Characteristics  = sam_d[sc]
        data_json.append({'Characteristics':Characteristics[0], 'Title':Title[0], 'Accesion_number':Accesion_number,
                          'x':str(data_2d[i][0]), 'y':str(data_2d[i][1]), 'similar': sim})
    else:
        data_json.append({'Characteristics': 'Unknown', 'Title': sc, 'Accesion_number': gse,'x': str(data_2d[i][0]),
                          'y':str(data_2d[i][1]), 'similar':sim})
    i+=1



with open('./data.json', 'w') as f:
    json.dump(data_json,f)
