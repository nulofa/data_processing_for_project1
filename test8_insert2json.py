import json
with open('./data.json', 'r') as f:
    jdata = json.load(f)
list = [1,3,4,5,6,7,9]
# with open('./data.json', 'w') as f:
#     for i in list:
#         print(jdata[i])