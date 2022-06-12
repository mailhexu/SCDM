import json

# test = dict(int_var=3, float_var=3.5, string_var="This is a string",
#            int_list_var=[1, 2, 3], float_list_var=[2.0, 3.0, 4.0])


test = dict(disentangle_func_type='unity', mu=0.1,
            sigma=0.5, exclude_bands=[1, 2, 3], nwann=3)

with open("input.json", 'w') as myfile:
    json.dump(test, myfile, indent=4)
