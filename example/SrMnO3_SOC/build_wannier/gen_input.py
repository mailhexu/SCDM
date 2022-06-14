import json

test = dict(method=1, disentangle_func_type='gauss', mu=-4,
            sigma=5, kmesh=[6, 6, 6],  nwann=5, project_to_anchor=True,
            anchor_kpt=[0, 0, 0], anchor_ibands=[46, 47, 48, 49, 50])

with open("input.json", 'w') as myfile:
    json.dump(test, myfile, indent=4)
