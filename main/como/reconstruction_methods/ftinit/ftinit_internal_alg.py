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
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, hstack, vstack
from scipy.sparse import eye as speye

####---------------------------------optimize_prob_scipy


def optimize_prob_scipy(prob, params, verbose=False):
    # SciPy implementation of the MILP solver

    # Convert to dense arrays if needed
    c = prob["c"].toarray().flatten() if hasattr(prob["c"], "toarray") else prob["c"]
    A = prob["A"].toarray() if hasattr(prob["A"], "toarray") else prob["A"]

    # Handle constraint senses (E=equality, L=<=, G>=)
    if isinstance(prob["csense"], str):
        # All constraints same type
        if prob["csense"] == "E":
            constraints = LinearConstraint(A, lb=prob["b"], ub=prob["b"])
        elif prob["csense"] == "L":
            constraints = LinearConstraint(A, lb=-np.inf, ub=prob["b"])
        else:  # 'G'
            constraints = LinearConstraint(A, lb=prob["b"], ub=np.inf)
    else:
        # Mixed constraints - process each type
        constraints = []
        for sense, row, b in zip(prob["csense"], A, prob["b"]):
            if sense == "E":
                constraints.append(LinearConstraint([row], lb=b, ub=b))
            elif sense == "L":
                constraints.append(LinearConstraint([row], lb=-np.inf, ub=b))
            else:  # 'G'
                constraints.append(LinearConstraint([row], lb=b, ub=np.inf))

        # Combine constraints
        if constraints:
            constraints = sum(constraints[1:], constraints[0])

    # Variable types (0:continuous, 1:integer)
    integrality = np.where(np.array(prob["vartype"]) == "C", 0, 1)

    # Bounds
    bounds = Bounds(prob["lb"], prob["ub"])

    # Solve
    result = milp(
        c=c,
        constraints=constraints,
        integrality=integrality,
        bounds=bounds,
        options={
            "disp": verbose,
            "presolve": True,
            "time_limit": params.get("TimeLimit", None),
            "node_limit": params.get("NodeLimit", None),
            "tol": params.get("intTol", 1e-7),
        },
    )

    # Format result to match expected structure
    return {"x": result.x, "fun": result.fun, "status": result.status, "message": result.message, "success": result.success}


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
    n_rxns_with_on_off = n_rxns

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

    force_on_lim = 0.1
    force_on_lim_ess = np.minimum(np.abs(fluxes) * 0.99, 0.1)

    vars_per_neg_rev = 3

    # Figure out the number of variables needed for metabolomics
    if met_data is not None:
        n_met_vars = 2 * met_data.shape[0] + 4 * n_met_neg_rev + n_met_neg_irrev
    else:
        n_met_vars = 0

    s_row = csr_matrix((n_mets, n_pos_irrev * 2 + n_pos_rev * 7 + n_neg_irrev * 2 + n_neg_rev * (1 + vars_per_neg_rev) + n_ess_rev * 6 + n_met_vars))
    s_eye = speye(n_rxns)
    n_y_block = n_pos_irrev + n_pos_rev + n_neg_irrev + n_neg_rev

    pi_rows = hstack(
        [
            s_eye[pos_irrev_rxns, :],
            speye(n_pos_irrev) * -force_on_lim,
            csr_matrix((n_pos_irrev, n_pos_rev + n_neg_irrev + n_neg_rev)),
            speye(n_pos_irrev) * -1,
            csr_matrix((n_pos_irrev, n_ess_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    pr_rows1 = hstack(
        [
            s_eye[pos_rev_rxns, :],
            csr_matrix((n_pos_rev, n_y_block + n_pos_irrev)),
            speye(n_pos_rev) * -1,
            speye(n_pos_rev),
            csr_matrix((n_pos_rev, n_pos_rev * 4 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    pr_rows2 = hstack(
        [
            csr_matrix((n_pos_rev, n_rxns + n_pos_irrev)),
            speye(n_pos_rev) * -force_on_lim,
            csr_matrix((n_pos_rev, n_neg_irrev + n_neg_rev + n_pos_irrev)),
            speye(n_pos_rev),
            speye(n_pos_rev),
            csr_matrix((n_pos_rev, n_pos_rev)),
            speye(n_pos_rev) * -1,
            csr_matrix((n_pos_rev, n_pos_rev * 2 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    pr_rows3 = hstack(
        [
            csr_matrix((n_pos_rev, n_rxns + n_y_block + n_pos_irrev)),
            speye(n_pos_rev),
            csr_matrix((n_pos_rev, n_pos_rev)),
            speye(n_pos_rev) * -100,
            csr_matrix((n_pos_rev, n_pos_rev)),
            speye(n_pos_rev),
            csr_matrix((n_pos_rev, n_pos_rev + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    pr_rows4 = hstack(
        [
            csr_matrix((n_pos_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev)),
            speye(n_pos_rev),
            speye(n_pos_rev) * 100,
            csr_matrix((n_pos_rev, n_pos_rev * 2)),
            speye(n_pos_rev),
            csr_matrix((n_pos_rev, n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    ni_rows = hstack(
        [
            s_eye[neg_irrev_rxns, :],
            csr_matrix((n_neg_irrev, n_pos_irrev + n_pos_rev)),
            speye(n_neg_irrev) * -100,
            csr_matrix((n_neg_irrev, n_neg_rev + n_pos_irrev + n_pos_rev * 6)),
            speye(n_neg_irrev),
            csr_matrix((n_neg_irrev, n_neg_rev * vars_per_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    nr_rows1 = hstack(
        [
            s_eye[neg_rev_rxns, :],
            csr_matrix((n_neg_rev, n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev)),
            speye(n_neg_rev) * -1,
            speye(n_neg_rev),
            csr_matrix((n_neg_rev, n_neg_rev + n_ess_rev * 6 + n_met_vars)),
        ]
    )

    nr_rows2 = hstack(
        [
            csr_matrix((n_neg_rev, n_rxns + n_pos_irrev + n_pos_rev + n_neg_irrev)),
            speye(n_neg_rev) * -100,
            csr_matrix((n_neg_rev, n_pos_irrev + n_pos_rev * 6 + n_neg_irrev)),
            speye(n_neg_rev),
            speye(n_neg_rev),
            speye(n_neg_rev),
            csr_matrix((n_neg_rev, n_ess_rev * 6 + n_met_vars)),
        ]
    )

    nr_rows = vstack([nr_rows1, nr_rows2])

    er_rows1 = hstack(
        [
            s_eye[ess_rev_rxns, :],
            csr_matrix((n_ess_rev, n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev)),
            speye(n_ess_rev) * -1,
            speye(n_ess_rev),
            csr_matrix((n_ess_rev, n_ess_rev * 4 + n_met_vars)),
        ]
    )

    er_rows2 = hstack(
        [
            csr_matrix((n_ess_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev)),
            speye(n_ess_rev),
            speye(n_ess_rev),
            speye(n_ess_rev) * -1,
            csr_matrix((n_ess_rev, n_ess_rev * 3 + n_met_vars)),
        ]
    )

    er_rows3 = hstack(
        [
            csr_matrix((n_ess_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev)),
            speye(n_ess_rev),
            csr_matrix((n_ess_rev, n_ess_rev * 2)),
            speye(n_ess_rev) * -100,
            speye(n_ess_rev),
            csr_matrix((n_ess_rev, n_ess_rev + n_met_vars)),
        ]
    )

    er_rows4 = hstack(
        [
            csr_matrix((n_ess_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev)),
            speye(n_ess_rev),
            csr_matrix((n_ess_rev, n_ess_rev)),
            speye(n_ess_rev) * 100,
            speye(n_ess_rev),
            csr_matrix((n_ess_rev, n_ess_rev)),
        ]
    )

    ##--------------------------------start of the mets--------------------------------##
    if met_data is not None:
        # Order of rxn "on" vars: Ypi Ypr Yni Ynr
        srt_met_data = csr_matrix(np.hstack([met_data[:, pos_irrev_rxns], met_data[:, neg_irrev_rxns], met_data[:, neg_rev_rxns]]))

        # First the setup for giving bonus if the met is included
        met_rows1 = hstack(
            [
                csr_matrix((n_metabol_mets, n_rxns)),
                -srt_met_data,
                csr_matrix((n_metabol_mets, n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6)),
                speye(n_metabol_mets),
                speye(n_metabol_mets),
                csr_matrix((n_metabol_mets, 4 * n_met_neg_rev + n_met_neg_irrev)),
            ]
        )

        # Then negative rev
        nr_eye = speye(n_neg_rev)

        met_rows2 = hstack(
            [
                csr_matrix((n_met_neg_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev)),
                nr_eye[met_neg_rev[neg_rev_rxns], :],
                csr_matrix((n_met_neg_rev, n_neg_rev * (vars_per_neg_rev - 1) + n_ess_rev * 6 + n_metabol_mets * 2)),
                speye(n_met_neg_rev) * -100,
                speye(n_met_neg_rev),
                csr_matrix((n_met_neg_rev, n_met_neg_rev * 2 + n_met_neg_irrev)),
            ]
        )

        met_rows3 = hstack(
            [
                csr_matrix((n_met_neg_rev, n_rxns + n_y_block + n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev)),
                nr_eye[met_neg_rev[neg_rev_rxns], :],
                csr_matrix((n_met_neg_rev, n_neg_rev * (vars_per_neg_rev - 2) + n_ess_rev * 6 + n_metabol_mets * 2)),
                speye(n_met_neg_rev) * 100,
                csr_matrix((n_met_neg_rev, n_met_neg_rev)),
                speye(n_met_neg_rev),
                csr_matrix((n_met_neg_rev, n_met_neg_rev + n_met_neg_irrev)),
            ]
        )

        met_rows4 = hstack(
            [
                csr_matrix((n_met_neg_rev, n_rxns + n_pos_irrev + n_pos_rev + n_neg_irrev)),
                nr_eye[met_neg_rev[neg_rev_rxns], :] * -0.1,
                csr_matrix((n_met_neg_rev, n_pos_irrev + n_pos_rev * 6 + n_neg_irrev)),
                nr_eye[met_neg_rev[neg_rev_rxns], :],
                nr_eye[met_neg_rev[neg_rev_rxns], :],
                csr_matrix((n_met_neg_rev, n_neg_rev * (vars_per_neg_rev - 2) + n_ess_rev * 6 + n_metabol_mets * 2 + n_met_neg_rev * 3)),
                speye(n_met_neg_rev) * -1,
                csr_matrix((n_met_neg_rev, n_met_neg_irrev)),
            ]
        )

        # Then negative irrev, i.e. flux -0.1*Yi -vnim == 0
        ni_eye = speye(n_neg_irrev)
        met_rows5 = hstack(
            [
                s_eye[met_neg_irrev, :],
                csr_matrix((n_met_neg_irrev, n_pos_irrev + n_pos_rev)),
                ni_eye[met_neg_irrev[neg_irrev_rxns], :] * -0.1,
                csr_matrix(
                    (
                        n_met_neg_irrev,
                        n_neg_rev
                        + n_pos_irrev
                        + n_pos_rev * 6
                        + n_neg_irrev
                        + n_neg_rev * vars_per_neg_rev
                        + n_ess_rev * 6
                        + n_metabol_mets * 2
                        + 4 * n_met_neg_rev,
                    )
                ),
                -speye(n_met_neg_irrev),
            ]
        )

        met_rows = vstack([met_rows1, met_rows2, met_rows3, met_rows4, met_rows5])
        met_var_c = np.concatenate([-np.ones(n_metabol_mets) * prod_weight, np.zeros(n_met_vars - n_metabol_mets)])

        # Variable bounds and types
        met_lb = np.concatenate(
            [np.zeros(n_metabol_mets * 2 + 2 * n_met_neg_rev), -100 * np.ones(n_met_neg_rev), np.zeros(n_met_neg_rev + n_met_neg_irrev)]
        )

        met_ub = np.concatenate(
            [np.ones(n_metabol_mets), np.inf * np.ones(n_metabol_mets), np.ones(n_met_neg_rev), np.inf * np.ones(3 * n_met_neg_rev + n_met_neg_irrev)]
        )

        met_var_type = np.concatenate(
            [np.array(["C"] * (n_metabol_mets * 2)), np.array(["B"] * n_met_neg_rev), np.array(["C"] * (3 * n_met_neg_rev + n_met_neg_irrev))]
        )

    else:
        met_rows = []
        met_var_c = []
        met_lb = []
        met_ub = []
        met_var_type = []

    # Prepare the optimization problem in a dictionary format
    prob = {}

    prob["a"] = vstack([s_row, pi_rows, pr_rows1, pr_rows2, pr_rows3, pr_rows4, ni_rows, nr_rows, er_rows1, er_rows2, er_rows3, er_rows4, met_rows])
    prob["b"] = np.zeros(prob["a"].shape[0])
    prob["c"] = np.concatenate(
        [
            np.zeros(n_rxns),
            -rxn_scores[pos_irrev_rxns],
            -rxn_scores[pos_rev_rxns],
            -rxn_scores[neg_irrev_rxns],
            -rxn_scores[neg_rev_rxns],
            np.zeros(n_pos_irrev + n_pos_rev * 6 + n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 6),
            met_var_c,
        ]
    )
    prob["lb"] = np.concatenate(
        [
            milp_model["lb"],
            np.zeros(n_y_block + n_pos_irrev + 5 * n_pos_rev),
            -100 * np.ones(n_pos_rev),
            np.zeros(n_neg_irrev + n_neg_rev * vars_per_neg_rev + n_ess_rev * 2),
            force_on_lim * np.ones(n_ess_rev),
            np.zeros(2 * n_ess_rev),
            -100 * np.ones(n_ess_rev),
            met_lb,
        ]
    )
    prob["lb"][ess_irrev_rxns] = force_on_lim_ess[ess_irrev_rxns]  # force flux for the essential irreversible reactions

    prob["ub"] = np.concatenate(
        [
            milp_model["ub"],
            np.ones(n_y_block),
            np.inf * np.ones(n_pos_irrev + n_pos_rev * 2),
            -100 * np.ones(n_pos_rev),
            np.inf * np.ones(3 * n_pos_rev + n_neg_irrev + vars_per_neg_rev * n_neg_rev + 3 * n_ess_rev),
            np.ones(n_ess_rev),
            np.inf * np.ones(2 * n_ess_rev),
            met_ub,
        ]
    )
    prob["vartype"] = np.concatenate(
        [
            np.array(["C"] * (n_rxns + n_pos_irrev + n_pos_rev)),
            np.array(["B"] * (n_neg_irrev + n_neg_rev)),
            np.array(["C"] * (n_pos_irrev + 2 * n_pos_rev)),
            np.array(["B"] * n_pos_rev),
            np.array(["C"] * (3 * n_pos_rev + n_neg_irrev + vars_per_neg_rev * n_neg_rev + 3 * n_ess_rev)),
            np.array(["B"] * n_ess_rev),
            np.array(["C"] * (2 * n_ess_rev)),
            met_var_type,
        ]
    )

    onoff_var_ind = np.arange(1, n_pos_irrev + n_pos_rev + n_neg_irrev + n_neg_rev + 1) + n_rxns
    onoff_pos_irrev = onoff_var_ind[:n_pos_irrev]
    onoff_pos_rev = onoff_var_ind[n_pos_irrev : n_pos_irrev + n_pos_rev]
    onoff_neg_irrev = onoff_var_ind[n_pos_irrev + n_pos_rev : n_pos_irrev + n_pos_rev + n_neg_irrev]
    onoff_neg_rev = onoff_var_ind[n_pos_irrev + n_pos_rev + n_neg_irrev :]
    met_var_ind = np.arange(1, n_metabol_mets + 1) + (len(prob["var_type"]) - n_met_vars)

    met_var_ind = np.arange(1, n_metabol_mets + 1) + (len(prob["var_type"]) - n_met_vars)
    if params.get("allowExcretion", False):
        prob["csense"] = np.concatenate(
            [np.array(["L"] * len(milp_model["mets"])), np.array(["E"] * (len(prob["b"]) - len(milp_model["mets"])))]
        ).tolist()
    else:
        prob["csense"] = "="
    params["intTol"] = 10**-7  # This value is very important. If set too low there is a risk that gurobi fails due to numerical issues -
    # this happened for Gurobi v. 10.0 with TestModelL. On the other hand, it shouldn't be too large either. With this value, fluxe fluxes of 10-7
    # can slip through, which should be fin, Another option if this becomes a problem is to set NumericFocus = 2, which makes the solver slower but fixes the issue with TestModelL.

    prob["osense"] = 1  # minimize
    if "startVals" in params and params["startVals"] is not None:
        prob["start"] = params["startVals"]  # This doesn't work...
    res = optimize_prob_scipy(prob, params, verbose)

    """
    Cannot find check solution function within ftinit
    # if not check_solution(res):
    #     if res["origStat"] == "TIME_LIMIT":
    #         em = "Time limit reached without finding a solution. Try increasing the TimeLimit parameter."
    #     else:
    #         em = "The problem is infeasible"
    #     print(em)  # Replace dispEM with print for simplicity
    """
    # Get the on/off vals
    onoff = np.ones(n_rxns_with_on_off)  # standard is on
    onoff[pos_irrev_rxns] = res["full"][onoff_pos_irrev]
    onoff[pos_rev_rxns] = res["full"][onoff_pos_rev]
    onoff[neg_irrev_rxns] = res["full"][onoff_neg_irrev]
    onoff[neg_rev_rxns] = res["full"][onoff_neg_rev]
    onoff2 = np.zeros(n_rxns_with_on_off)  # standard is off
    onoff2[pos_irrev_rxns] = res["full"][onoff_pos_irrev]
    onoff2[pos_rev_rxns] = res["full"][onoff_pos_rev]
    onoff2[neg_irrev_rxns] = res["full"][onoff_neg_irrev]
    onoff2[neg_rev_rxns] = res["full"][onoff_neg_rev]
    # Get all reactions used in the irreversible model
    deleted_rxns = (onoff < 0.5) & (rxn_scores != 0)
    turned_on_rxns = (onoff2 >= 0.5) & (rxn_scores != 0)
    fluxes = res["full"][:n_rxns]
    # Extract the met data
    met_production = np.round(res["full"][met_var_ind]).astype(bool)

    return deleted_rxns, met_production, res, turned_on_rxns, fluxes
