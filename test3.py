import numpy as np

sam_d = np.load('./sam_details_new1.npy', allow_pickle=True)
sam_d=sam_d.item()
sam_d.pop('GSE77211fold_change(ago/wt)')
np.save('./sam_details_new1.npy', sam_d)
# sam = np.load('./sam_details2.npy', allow_pickle=True)
# sam = sam.item()
# for s in sam:
#     print(s, sam[s])