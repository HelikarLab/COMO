# This is the Python version of ftINITInternalAlg.m (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITInternalAlg.m

# ftINITInternalAlg: This function runs the MILP for a step in ftINIT

## INPUT
# model: a reference model structure
# rxn_scores: a vector of scores for the reactions in the model.
#             positive scores are reactions to keep and negative scores are reactions to exclude. rxns set to 0 are excluded from the problem.
# met_data: boolean matrix with mets as rows and rxns as columns saying which reaction produces each detected met (optional, default [])
# essential_rxns: cell array of reactions that are essential and that have to be in the resulting model. this is normally
#                 used when fitting a model to task (fee fitTasks)
# prod_weight: a score that determines the value of having net-production of metabolites. this is a way of having a more functional
#              network as it provides a reason on for including bad reactions for connectivity reasons. this score is for each metabolite,
#              and the sum of these weights and the scores for the reactions is what is optimized
# allow_excretion: true if excretion of all metabolites should be allowed. this results in fewer reactions being considered dead-ends,
#                  but all reactions in the resulting model may not be able to carry flux. if this is "false" then the equality constraints
#                  are taken from model.b. if the input model lacks exchange reactions then this should probably be "true", or a large proportion
#                  of the model would be excluded for connectivity reasons.
# rem_posrev: if true, the positive reversible reactions are removed from the problem. this is used in step 1 of ftINIT.
# params: parameters for the MILP, for example MIPGap and TimeLimit
# start_vals: Start values for the MILP, typically used when rerunning with a higher MIPGap, to use the results from the previous run
# fluxes: fluxes from the last run
# verbose: if true, the MILP progression will be shown.

## OUTPUT
# deleted_rxns: reactions which were deleted by the algorithm (only rxns included in the problem)
# met_production: array that indicates which of the metabolites in presentMets that could be produced
#                 0: metabolite could not be produced
#                 1: metabolite could be produced
# res: the result from the MILP
# turned_on_rxns: the reactions determined to be present (only rrxns included in the problem)
# fluxes: the fluxes from the MILP

# This function is the actual implementation of the algorithm. See ftINIT for a higher-level function for model reconstruction.

import numpy as np


def ftinit_internal_alg(model, rxn_scores, met_data, essential_rxns, prod_weight, allow_excretion, rem_pos_rev, params, start_vals, fluxes, verbose):
    essential_rxns = params.get("essential_rxns", [])
    essential_rxns = np.array(essential_rxns).flatten()
    prod_weight = params.get("prod_weight", 0.5)

    # The model should be in the reversible format and all relevant exchange reactions should be open
    if "constrained" in model:
        print("Exchange metabolites are still present in the model. Use simplifyModel if this is not intended")
    essential = np.isin(model.reactions, essential_rxns)

    # Some nice to have numbers
    n_mets = len(model.metabolites)
    n_rxns = len(model.reactions)
    nrxns_with_on_off = n_rxns

    # Reactions with score 0 will just be left in the model, and will not be a part of the problem (but can carry flux).
    # It is possible to set score = 0 for e.g. spontaneous reactions, exchange rxns, etc., which may not be that interesting to remove

    # If makeIrrev is on, we can just as well skip all positive reversible rxns - they will be turned on without carrying
    # any flux since they can form a loop within themself (fwd-rev)
    rem_pos_rev = params.get("rem_pos_rev", False)
    if rem_pos_rev:
        rxn_scores[(rxn_scores > 0) & (model["rev"] != 0)] = 0

    # Handle metabolomics
    # A problem with the metabolomics is that some of the producer reactions for a metabolite could be excluded from the problem.
    # We solve this by adding them with score 0. They can still be seen as positive, since they either don't matter (there is another producer)
    # or they can contribute to an increased score.
    rev_rxns = model.rev != 0
    ess_rev_rxns = np.where(rev_rxns & essential)[0]
    ess_irrev_rxns = np.where(~rev_rxns & essential)[0]

    met_data = params.get("met_data", None)
    if met_data is not None:
        # Remove any metData rows that are connected to an essential rxn
        # these will be on regardless and will cause prolems below
        contains_essential = np.any(met_data[:, ess_rev_rxns], axis=1)
        met_data = met_data[~contains_essential, :]

    if met_data is not None:
        met_rxns = np.any(met_data, axis=0)
        pos_rxns = (rxn_scores > 0) | ((rxn_scores == 0) & met_rxns)
    else:
        pos_rxns = rxn_scores > 0

    neg_rxns = rxn_scores < 0

    pos_rev_rxns = np.where(pos_rxns & rev_rxns & ~essential)[0]
    neg_rev_rxns = np.where(neg_rxns & rev_rxns & ~essential)[0]
    pos_irrev_rxns = np.where(pos_rxns & ~rev_rxns & ~essential)[0]
    neg_irrev_rxns = np.where(neg_rxns & ~rev_rxns & ~essential)[0]

    n_pos_rev = len(pos_rev_rxns)
    n_neg_rev = len(neg_rev_rxns)
    n_pos_irrev = len(pos_irrev_rxns)
    n_neg_irrev = len(neg_irrev_rxns)
    n_ess_rev = len(ess_rev_rxns)
    n_ess_irrev = len(ess_irrev_rxns)  # not used byt left for symmetry
    n_metabol_mets = met_data.shape[0] if met_data is not None else 0

    if met_data is not None:
        met_neg_rev = neg_rxns & rev_rxns & ~essential & met_rxns
        met_neg_irrev = neg_rxns & ~rev_rxns & ~essential & met_rxns
        n_met_neg_rev = np.sum(met_neg_rev)
        n_met_neg_irrev = np.sum(met_neg_irrev)

    milp_model = model

    # These six categories need to be handled separately:
    #
    # pos_irrev (positive reaction score, irreversible):
    # flux >= 0.1*Yi, flux - 0.1*Yi - VPI == 0, 0 <= Yi(onoff) <= 1, 0 <= VPI <= Inf
    # The nice thing with rxns with positive score is that they do not require a boolean. The optimizer will strive for maximizing
    # the Yi here, so it can be a continuous variable, where we just force a flux it the variable is on.
    #
    # pos_rev:
    # 1. split up the flux into positive and negative: flux - vnrp + vnrn == 0, 0 <= vprp,vprn <= Inf.
    # 2. force one of the fluxes on: vprp + vprn >= 0.1*Yi, vprp + vprn - 0.1*Yi -vprv1 == 0. 0 <= Yi(onoff) <= 1, vprv1 >= 0
    # 3. we need a bool to make sure that one of vprp and vprn are always zero:
    #    vprb E {0,1}: vprp <= 100*vprb (if bool is 0, vprp is zero): Vprp - 100*vprb + vprv2 == 0, vprv2 >= 0
    #    vprn <= (1-vprb)*100 (if bool is one, vprn is zero): vprn + 100*vprb + vprv3 == 0, -100 <= vprv3 < inf
    #
    # neg_irrev:
    # for negative scores, we need to use an integer, I see no way around that.
    # flux < 100*Yi, Yi E {0,1}. flux - 100*Yi + vni == 0
    #
    # neg_rev:
    # two cases:
    # A: without irrev model
    #   1: split up the flux into positive and negative: flux - vnrp + vnrn == 0, 0 <= vprp,vprn <= Inf.
    #   2: Force the Yi (on/off) var on if the flux is on: vnrp + vnrn <= 0.1*Yi, Yi E {0,1}: vnrp + vnrn - 0.1 * Yi + vnrv1 == 0, vnrv1 >= 0.
    # B: irrev model (Not used for now)
    #   1: we now have two reactions, but still just boolean variable. we know that the flux can onl be positive, so we just say that
    #   flux1 + flux2 <= 0.1 * Yi, Yi E {0,1}: flux1 + flux2 - 0.1 * Yi + vnrv1 == 0, vnrv >= 0. we don't care if there is any loop - there is just
    #   no benefit for the objective to generate one, and it doesn't matter.
    #
    # ess_irrev (essential rxn, i.e. forced to carry flux, irreversible):
    # flux >= 0.1 - solved with lb = 0.1
    #
    # ess_rev (these are not really used):
    # 1. split up the flux intp positive and negative: flux - verp + vern == 0, 0 <= verp,vern <= Inf
    # 2. force one of the fluxes on: verp + vern >= 0.1, verp + vern - verv1 == 0. verv1 >= 0.1
    # 3. we need a bool to make sure that one of verp and vern are always zero:
    #   verb E {0,1}: verp <= 100*verb (if bool is 0, verp is zero): verp - 100*verb + verv2 == 0, verv2 >0
    #   vern <= (1-verb)*100 (if bool is one, vern is zero): vern + 100*verb + verv3 == 0, -100 <= verv3 <= inf

    # Now metabolomics
    # the basic idea is that the bonus is gotten if any of the producer reactions are on. we can therefore use the "one" variable directly for
    # the reactions that are included. we therefore add one variable per met, mon, which can be continuous and between 0 and 1. we then can say that
    # mon <= v1 + v2 + ... + vn, i.e. that mon + mv1 - v1 - v2 ... - vn == 0, 0 <= mon <= 1, 0 <= mv1 <= inf
    # mon gets -prodweight in the c vector (i.e. objective function_
    # a tricky part is that some of the reactions producing a metabolite are left outside the problem. this is solved by moving them into the
    # problem, but with the score 0. they can still be treated as positive reactions. in the end, we are not interested if they are on though, but
    # they may allow for production of a metabolite, giving no benefit of turning on another producer reaction. another tricky thing is that some
    # metabolite producers have negative score. there is no automatic mechanism for forcing fluc on if the variable is on - this is simply not needed
    # for negative variables -  unless they are a metabolite producer. we therefore need to add an extra constraint there, similar to the case of
    # positive. it gets a bit complicated.
    # for reversible, we need a bool to make sure that one of vnrp and vnrn stays zero:
    #   vnrbm E {0,1}: vnrp <= 100*vnrbm (if bool is 0, vnrp is zero): vnrp - 100*vnrbm + vnrvm1 == 0, vnrvm1 >= 0
    #   vnrn <= (1-vnrbm)*100 (if bool is one, vnrn is zero): vnrn + 100*vnrbm + vnrvm2 == 0, -100 <= vnrvm2 <= inf
    # we then also say that vnrp + vnrn >= 0.1*Yi, vnrp + vnrn - 0.1*Yi - vnrvm3 == 0, vnrvm3 >= 0
    # for irreversible, it is fairly straight-forward. we just day that flux >= 0.1*Yi, i.e. flux - 0.1*Yi - vnim == 0, 0<= vnim <= Inf

    # the total matrix then looks like this. note that we have ordered the categories, two reactions each, in the S matrix to make it easier to
    # follow the figure. in practice, they come in random order, which is why the SEye variable is used further down.
    # all categories below have two columns - the two column starts where label starts
    # since the label cannot fit in two chars, they are on multiple lines
    #
    # the met variables and constraints are not in this figure for practical reasons
    #
    #                           vprb  vprv3 vnrn  vern  verv2 mv1
    # pi  ni  ei  Ypi Yni vpi vprn  vprv2 vnrp  verp  verb  mon
    #   pr  nr  er  Ypr Ynr vprp  vprv1 vni   vnrv1 verv1 verv3
    # SSSSSSSSSSSS0000000000000000000000000000000000000000000000
    # SSSSSSSSSSSS0000000000000000000000000000000000000000000000 S block
    # SSSSSSSSSSSS0000000000000000000000000000000000000000000000
    # 1           N       -                                      Pos irrev block
    #  1           N       -
    #   1                   - 1                                  Pos rev block 1
    #    1                   - 1
    #               N       1 1   -                              Pos rev block 2
    #                N       1 1   -
    #                       1   M   1                            Pos rev block 3
    #                        1   M   1
    #                         1 C     1                          Pos rev block 4
    #                          1 C     1
    #     1           M                 1                        Neg irrev block
    #      1           M                 1
    #       1                             - 1                    Neg rev block 1
    #        1                             - 1
    #                   M                 1 1 1                  Neg rev block 2
    #                    M                 1 1 1
    #           1                               - 1              Ess rev block 1
    #            1                               - 1
    #                                           1 1 -            Ess rev block 2
    #                                            1 1 -
    #                                           1     M 1        Ess rev block 3
    #                                            1     M 1
    #                                             1   C   1      Ess rev block 4
    #                                              1   C   1
    #             --------                                  1 1  Met block - Here, we assume that all variables support each met.
    #             --------                                   1 1             In practice, fewer of the "on" variables with -1 should be included
    #
    #
    # M = -100
    # N = -0.1
    # D = 0.1
    # C = 100
    # - = -1
    #
    # It looks a bit different for the irrev model, but still quite similar
    # build the A matrix
    # S row
    # when forcing on essential rxns, use the flux value of the previous run (set to 0.1 the first time)
    # don't set it above 0.1, may starve something else out. leave a margin of 1% from  the last run.

    return deleted_rxns, met_production, res, turned_on_rxns, fluxes
