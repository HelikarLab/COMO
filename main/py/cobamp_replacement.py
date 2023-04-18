import cobra
import numpy as np


# from cobamp.wrappers.external_wrappers import COBRAModelObjectReader

def get_reaction_and_metabolite_ids(cobra_model: cobra.Model) -> tuple:
    return tuple(
        [
            [x.id for x in lst] for lst in (cobra_model.reactions, cobra_model.metabolites)
        ]
    )


def get_stoichiometric_matrix(cobra_model: cobra.Model, reaction_ids: list, metabolite_ids: list) -> np.array:
    """
    This function will retrieve the stoichiometric matrix from a cobra model.

    This is used instead of cobamp's get_stoichiometric_matrix() function to avoid dependencies on cobamp
    Inspiration taken from: https://github.com/BioSystemsUM/cobamp/blob/master/src/cobamp/wrappers/cobra.py

    """
    matrix = np.zeros((len(metabolite_ids), len(reaction_ids)))
    for i, reaction_id in enumerate(reaction_ids):
        for metabolite, coefficient in cobra_model.reactions.get_by_id(reaction_id).metabolites.items():
            matrix[metabolite_ids.index(metabolite.id), i] = coefficient
    return matrix


def get_model_bounds(cobra_model: cobra.Model) -> tuple[list[float], list[float]]:
    model_bounds = [reaction.bounds for reaction in cobra_model.reactions]
    lower_bounds = []
    upper_bounds = []
    for lower_bound, upper_bound in model_bounds:
        lower_bounds.append(lower_bound)
        upper_bounds.append(upper_bound)
    return lower_bounds, upper_bounds


def main():
    model_file: str = "/Users/joshl/PycharmProjects/MADRID/main/data/GeneralModelUpdatedV2.json"
    cobra_model: cobra.Model = cobra.io.load_json_model(model_file)
    
    cobra_model.solver = "gurobi"
    reaction_ids, metabolite_ids = get_reaction_and_metabolite_ids(cobra_model)
    stoichiometric_matrix: np.array = get_stoichiometric_matrix(cobra_model, reaction_ids, metabolite_ids)
    lower_bounds, upper_bounds = get_model_bounds(cobra_model)
    
    print("DONE")


if __name__ == '__main__':
    main()
