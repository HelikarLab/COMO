from dataclasses import dataclass
from typing import Annotated, Literal

import cobra
import numpy as np
import scipy


@dataclass(frozen=True, kw_only=True, slots=True)
class FisherExactTest:
    """Evalute Fisher's Exact Test for reaction presence.

    Fisher's Exact Test is a non-parametric statistical test used to determine if there are nonrandom associations between two variables.
    It is useful in metabolic modeling because it can help assess whether the presence or absence of certain reactions in a metabolic model
        is independent of a specific condition or treatment without assuming the distribution of the data.

    To calculate the Fisher's Exact Test, execute :func:`FisherExactTest.run`, which will return an instance of :class:`FisherExactTest`

    References:
        [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)
        [Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_exact_test)
    """

    pathway: Annotated[str, "The pathway test"]
    statistic: Annotated[float, "The odds ratio of the test"]
    pvalue: Annotated[float, "The p-value of the test"]
    a: Annotated[int, "Number of reactions in the pathway AND scenario model"]
    b: Annotated[int, "Number of reactions in the pathway but NOT the scenario model"]
    c: Annotated[int, "Number of reactions NOT in the pathway but ARE in the scenario model"]
    d: Annotated[int, "Number of reactions NOT in the pathway NOR the scenario model"]

    @classmethod
    def run(
        cls: type["FisherExactTest"],
        reference: cobra.Model,
        scenario: cobra.Model,
        pathway: str,
        alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    ) -> "FisherExactTest":
        """Perform a Fisher's Exact Test on two models with a known reference model.

        This test is based on the following assumptions:
        - The general "reference" model was used to reconstruct the comprehensive and scenario models (such as Recon3D)
        - A scenario-specific model exists that may not be representative of true biology

        ---

        Given the following contingency table for a set of conditions and N reactions:
        - A: Reactions in `pathway` and the scenario-specific model
        - B: Reactions in `pathway` but not the scenario-specific model
        - C: Reactions in the scenario-specific model that are not a part of `pathway`
        - D: Reactions not in `pathway` that are also not found in the scenario-specific model

        | Reaction Status  | In scenario-specific model | Not in scenario-specific model |          Row Total           |
        |:----------------:|:--------------------------:|:------------------------------:|:----------------------------:|
        |   In `pathway`   |             A              |               B                |            A + B             |
        | Not in `pathway` |             C              |               D                |            C + D             |
        |   Column Total   |           A + C            |             B + D              | A + B + C + D (=N reactions) |

        A two-sided Fisher's exact test will ask the question:
            > Is the inclusion or exclusion of this reaction in the patient model independent of its status in the reference model?

        If the scenario-specific dataset is "small", the reconstruction will likely have excluded/dropped some reactions
            as a result of the limited data available. This means the Fisher's exact test may show **many apparent differences**.
            However, this could be noise from undersampling and not indicative of the true underlying biology.
        In practice, if only a few reactions fall into "condition A" (above), this suggests that the scenario-specific model is too sparse
            and not reconstructed with enough data.

        Args:
                reference: The general reference model that was used to build the model (e.g., Recon3D)
                scenario: The scenario-specific model to test (e.g., built using a small cohort of single-cell RNA-seq data)
                pathway: The pathway to investigate for a Fisher's Exact Test
                alternative: The alternative hypothesis to test

        Returns:
                The p-value indicating whether the reaction presence in the scenario model is independent of the reference model.
        """
        scenario_rxn_ids: set[str] = {rxn.id for rxn in scenario.reactions}

        a = 0  # a reaction is in the given pathway and scenario model
        b = 0  # a reaction is in the given pathway but not the scenario model
        c = 0  # a reaction is not in the given pathway but is in the scenario model
        d = 0  # a reaction is not in the given pathway OR the scenario model

        for rxn in reference.reactions:
            in_pathway = rxn.subsystem == pathway
            in_scenario = rxn.id in scenario_rxn_ids

            if in_pathway:
                if in_scenario:
                    a += 1
                else:
                    b += 1
            else:
                if in_scenario:
                    c += 1
                else:
                    d += 1

        result = scipy.stats.fisher_exact(np.array([[a, b], [c, d]]), alternative=alternative)
        return cls(statistic=result.statistic, pvalue=result.pvalue, pathway=pathway, a=a, b=b, c=c, d=d)
