
dict = {"rxn1" : -1.,
        "rxn2" : -1.,
        "rxn3" :  0.,
        "rxn4" :  0.,
        "rxn5" :  1.,
        "rxn6" :  1.
        }

dict_exp = list({k:v for (k,v) in dict.items() if v > 0}.keys())



print(dict_exp)
